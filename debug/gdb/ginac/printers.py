
import ctypes
import fractions
import gdb
import gdb.printing
import gdb.types
import re
import sys


if sys.version_info.major == 3:
    long = int
    Iterator = object
else:
    class Iterator:
        def next(self):
            return self.__next__()


def cast_and_dereference(ptr):
    if ptr.type.code == gdb.TYPE_CODE_PTR:
        if ptr.type != ptr.dynamic_type:
            ptr = ptr.cast(ptr.dynamic_type)
        return ptr.dereference()


def unwrap_ex(val):
    basic_ptr = val['bp']['p']
    ptr = basic_ptr.cast(basic_ptr.dynamic_type)
    return ptr.dereference()


def ex_to_number(wrapped_val):
    val = unwrap_ex(wrapped_val)
    type_tag = val.type.unqualified().tag
    if type_tag == 'GiNaC::numeric':
        return val['value']


def as_a(ex, typename):
    val = unwrap_ex(ex)
    type_tag = val.type.unqualified().tag
    if type_tag == 'GiNaC::{}'.format(typename):
        return val


def find_type(orig, name):
    """Starting with orig, search for the member type NAME"""
    typ = orig.strip_typedefs()
    while True:
        search = '{0}::{1}'.format(typ.name, name)
        try:
            return gdb.lookup_type(search)
        except RuntimeError:
            pass
        # try the superclass
        field = typ.fields()[0]
        if not field.is_base_class:
            raise ValueError('Cannot find type %s::%s' % (str(orig), name))
        typ = field.type


class ExPrinter:

    def __init__(self, val):
        self.val = None
        ptr = val['bp']['p']
        if ptr.type.code == gdb.TYPE_CODE_PTR:
            dyntype = ptr.dynamic_type
            if dyntype != ptr.type:
                ptr = ptr.cast(dyntype)
            self.val = ptr.dereference()

    def to_string(self):
        return str(self.val)

    def display_hint(self):
        return 'ex'


class ClNumber(object):

    IMMEDIATE_MASK = 0b111
    VALUE_SHIFT = 3
    DIGIT_BITS = 64

    def __init__(self, val):
        self.val = val
        self.cl_class_bignum = gdb.parse_and_eval('::cln::cl_class_bignum')
        self.cl_class_ratio = gdb.parse_and_eval('::cln::cl_class_ratio')
        self.cl_class_complex = gdb.parse_and_eval('::cln::cl_class_complex')

    @property
    def word(self):
        return self.val['word']

    @property
    def pointer(self):
        return self.val['pointer']

    @property
    def is_immediate(self):
        return self.word & self.IMMEDIATE_MASK != 0

    @property
    def is_pointer(self):
        return not self.is_immediate

    @property
    def is_zero(self):
        return self.fixnum == 0

    def _decode_fixnum(self):
            return ctypes.c_long(self.word).value >> self.VALUE_SHIFT

    @property
    def fixnum(self):
        if self.is_immediate:
            return self._decode_fixnum()

    @property
    def is_bignum(self):
        heappointer = self.val['heappointer']
        heappointer = heappointer.cast(heappointer.dynamic_type)
        return heappointer.dereference()['type'] == self.cl_class_bignum.address

    @property
    def is_ratio(self):
        ptr = self.val['heappointer']
        return ptr.dereference()['type'] == self.cl_class_ratio.address

    @property
    def is_complex(self):
        ptr = self.val['heappointer']
        return ptr.dereference()['type'] == self.cl_class_complex.address

    # XXX: python complex are always floating-point
    @property
    def pynumber(self):
        if self.is_immediate:
            return self._decode_fixnum()
        elif self.is_bignum:
            return self.decode_bignum()
        elif self.is_ratio:
            return self.decode_ratio()

    @property
    def number(self):
        if self.is_immediate:
            return self._decode_fixnum()
        elif self.is_bignum:
            return self.decode_bignum()
        elif self.is_ratio:
            return self.decode_ratio()
        elif self.is_complex:
            return self.decode_complex()

    def decode_bignum(self):
        heappointer = self.val['heappointer']
        heappointer = heappointer.cast(heappointer.dynamic_type)
        heap_bignum_type = gdb.lookup_type('::cln::cl_heap_bignum')
        heap_bignum_pointer = heappointer.cast(heap_bignum_type.pointer())
        heap_bignum = heap_bignum_pointer.dereference()
        length = heap_bignum['length']
        unsigned_long_type = gdb.lookup_type('unsigned long')
        data_ptr = heap_bignum['data'].cast(unsigned_long_type.pointer())
        value = 0
        i = 0
        while i < length:
            value += int(data_ptr.dereference()) * 2**(i*self.DIGIT_BITS)
            data_ptr = data_ptr + 1
            i += 1
        return value

    def decode_ratio(self):
        ptr = self.val['heappointer']
        ratio_type = gdb.lookup_type('::cln::cl_heap_ratio')
        ratio_ptr = ptr.cast(ratio_type.pointer())
        ratio = ratio_ptr.dereference()
        numerator = ClNumber(ratio['numerator'])
        denominator = ClNumber(ratio['denominator'])
        return fractions.Fraction(numerator.number, denominator.number)

    def decode_complex(self):
        ptr = self.val['heappointer']
        complex_type = gdb.lookup_type('::cln::cl_heap_complex')
        complex_ptr = ptr.cast(complex_type.pointer())
        complex_ = complex_ptr.dereference()
        real = ClNumber(complex_['realpart']).number
        imag = ClNumber(complex_['imagpart']).number
        mark = '+I'
        if imag == 1:
            return str(real) + ' + I'
        elif imag == -1:
            return str(real) + ' - I'
        elif imag == 0:
            return str(real)
        elif imag < 0:
            return '{0} - {1}*I'.format(str(real), str(-imag))
        elif imag > 0:
            return '{0} + {1}*I'.format(str(real), str(imag))
        else:
            return '{0} + ({1})*I'.format(str(real), str(imag))


class ClNumberPrinter:

    def __init__(self, val):
        self.val = ClNumber(val)

    def to_string(self):
        return str(self.val.number)

    def display_hint(self):
        return 'number'


class NumericPrinter:

    def __init__(self, val):
        self.numval = val['value']

    def to_string(self):
        return ClNumberPrinter(self.numval).to_string()

    def display_hint(self):
        return 'numeric'


class ModIntegerPrinter:

    def __init__(self, val):
        self.val = val
        self.ring_type = gdb.lookup_type('cln::cl_heap_modint_ring')

    def to_string(self):
        ring_ptr = self.val['_ring']['pointer']
        ring = ring_ptr.cast(self.ring_type.pointer()).dereference()
        value = self.val['rep']
        modulus = ring['modulus']
        return '{0} mod {1}'.format(str(value), str(modulus))


def print_power(basis, exponent):
    fmt = '({basis})^({exponent})'
    basis_need_parens = True
    exponent_need_parens = True
    numeric_exponent = ex_to_number(exponent)
    if numeric_exponent:
        as_cl_N = ClNumber(numeric_exponent)
        as_fixnum = as_cl_N.fixnum
        if as_fixnum:
            exponent = as_fixnum
            if as_fixnum > 0:
                exponent_need_parens = False
    symbol_basis = as_a(basis, 'symbol')
    if symbol_basis:
        basis_need_parens = False
        basis = symbol_basis
    if basis_need_parens:
        fmt ='({basis})'
    else:
        fmt = '{basis}'
    if exponent == 1:
        return fmt.format(basis=str(basis))
    fmt = fmt + '^'
    if exponent_need_parens:
        fmt = fmt + '({exponent})'
    else:
        fmt = fmt + '{exponent}'
    return fmt.format(basis=str(basis), exponent=str(exponent))


class PowerPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        basis = self.val['basis']
        exponent = self.val['exponent']
        return print_power(basis, exponent)

    def display_hint(self):
        return 'power'


class StdVectorIterator(Iterator):

    def __init__(self, val):
        self.item = val['_M_impl']['_M_start']
        self.finish = val['_M_impl']['_M_finish']

    def __iter__(self):
        return self

    def __next__(self):
        if self.item == self.finish:
            raise StopIteration
        elt = self.item.dereference()
        self.item = self.item + 1
        return elt


def print_term(e, n = 0, leading_plus=True):
    term = e['rest']
    coeff = ex_to_number(e['coeff'])
    coeff_pynum = ClNumber(coeff).pynumber
    
    if coeff_pynum == 1:
        if n == 0 and not leading_plus:
            return str(term)
        else:
            return '+' + str(term)
    elif coeff_pynum == -1:
        return '-' + str(term)
    elif coeff_pynum > 0:
        if n == 0 and not leading_plus:
            fmt = '{0}*{1}'
        else:
            fmt = '+{0}*{1}'
    else:
        fmt = '{0}*{1}'
    return fmt.format(ClNumberPrinter(coeff).to_string(), str(term))


def print_sum(seq, leading_plus=False):
    return ''.join(print_term(elt, n, leading_plus=leading_plus)
                   for n, elt in enumerate(StdVectorIterator(seq)))

class SumPrinter:

    def __init__(self, val):
        self.val = val

    def should_skip_overall_coeff(self, oc):
        return ClNumber(oc['value']).is_zero

    def to_string(self):
        seq = self.val['seq']
        overall_coeff = unwrap_ex(self.val['overall_coeff'])
        if self.should_skip_overall_coeff(overall_coeff):
            return print_sum(seq)
        else:
            return str(overall_coeff) + print_sum(seq, leading_plus=True)

    def display_hint(self):
        return 'sum'


class ProdPrinter:

    def __init__(self, val):
        self.seq = val['seq']
        self.overall_coeff = ex_to_number(val['overall_coeff'])

    def should_skip_overall_coeff(self):
        return ClNumber(self.overall_coeff).fixnum == 1

    def print_overall_coeff(self):
        oc_as_pynum = ClNumber(self.overall_coeff).pynumber
        if oc_as_pynum:
            if oc_as_pynum == 1:
                return ''
            elif oc_as_pynum == -1:
                return '-'
            else:
                return str(oc_as_pynum) + '*'
        else:
            return str(self.overall_coeff) + '*'

    def print_term(self, term):
        base = term['rest']
        exponent = term['coeff']
        return print_power(base, exponent)

    def print_seq(self):
        return '*'.join([self.print_term(term) for term in
                         StdVectorIterator(self.seq)])

    def to_string(self):
        return self.print_overall_coeff() + self.print_seq()


class SymbolPrinter:
    "Print GiNaC::symbol"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        return str(self.val['name']).strip('"')

    def display_hint(self):
        return 'symbol'


class RelationalPrinter:

    def __init__(self, val):
        self.lhs = val['lh']
        self.rhs = val['rh']
        self.operator = val['o']

    def to_string(self):
        return '{0} {1} {2}'.format(str(self.lhs),
                                    str(self.operator),
                                    str(self.rhs))

    def display_hint(self):
        return 'relation'


class StdListIterator(Iterator):
    def __init__(self, nodetype, head):
        self.nodetype = nodetype
        self.head = head.address
        self.base = head['_M_next']
        self.count = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.base == self.head:
            raise StopIteration
        elt = self.base.cast(self.nodetype).dereference()
        self.base = elt['_M_next']
        count = self.count
        self.count = self.count + 1
        return ('[%d]' % count, elt['_M_data'])


class LstPrinter:

    def __init__(self, val):
        self.val = val['seq']

    def children(self):
        nodetype = find_type(self.val.type, '_Node')
        nodetype = nodetype.strip_typedefs().pointer()
        return StdListIterator(nodetype, self.val['_M_impl']['_M_node'])

    @property
    def is_empty(self):
        head = self.val['_M_impl']['_M_node']
        child = self.val['_M_impl']['_M_node']['_M_next']
        return head.address == child

    def to_string(self):
        if self.is_empty:
            return 'empty lst'
        else:
            return 'lst'


def str_lookup_function(val):
    lookup_tag = val.type.tag
    if val.type.code == gdb.TYPE_CODE_PTR:
        dyntype = val.dynamic_type
        if dyntype != val.type:
            val = val.cast(dyntype)
        lookup_tag = val.type.target().unqualified().tag
        val = val.dereference()
    if val.type.code == gdb.TYPE_CODE_REF:
        lookup_tag = val.type.target().unqualified().tag
        val = val.referenced_value()
    if lookup_tag is None:
        return None
    printers_table = {
        'GiNaC::ex': ExPrinter,
        'GiNaC::power': PowerPrinter,
        'GiNaC::add': SumPrinter,
        'GiNaC::mul': ProdPrinter,
        'GiNaC::numeric': NumericPrinter,
        'GiNaC::relational': RelationalPrinter,
        'GiNaC::container<std::__cxx11::list>': LstPrinter,
        'cln::cl_number': ClNumberPrinter,
        'cln::cl_N': ClNumberPrinter,
        'cln::cl_R': ClNumberPrinter,
        'cln::cl_RA': ClNumberPrinter,
        'cln::cl_I': ClNumberPrinter,
        'cln::cl_MI': ModIntegerPrinter,
    }
    printer = printers_table.get(lookup_tag)
    if printer:
        return printer(val)
    regex = re.compile('GiNaC::.*symbol')
    if regex.match(lookup_tag):
        return SymbolPrinter(val)
    return None


def register_libginac_printers(objfile):
    gdb.printing.register_pretty_printer(objfile, str_lookup_function)
