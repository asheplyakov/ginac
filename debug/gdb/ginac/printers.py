
import ctypes
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

    def __init__(self, val):
        self.val = val

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
        return self.word & self.IMMEDIATE_MASK == 0

    @property
    def is_zero(self):
        return self.fixnum == 0

    @property
    def fixnum(self):
        if self.is_immediate:
            return ctypes.c_long(self.word).value >> self.VALUE_SHIFT


class ClNumberPrinter:

    def __init__(self, val):
        self.val = ClNumber(val)

    def to_string(self):
        fn = self.val.fixnum
        if fn is not None:
            return str(fn)
        else:
            return str(self.val.pointer)

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


class PowerPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        return '(%s)^(%s)' % (str(self.val['basis']),
                              str(self.val['exponent']))

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


def print_term(e):
    term = e['rest']
    coeff = ex_to_number(e['coeff'])
    if ClNumber(coeff).fixnum == 1:
        return str(term)
    else:
        return '{0}*{1}'.format(ClNumberPrinter(coeff).to_string(), str(term))


def print_sum(seq):
    return ' + '.join([print_term(elt) for elt in StdVectorIterator(seq)])


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
            return str(overall_coeff) + ' + ' + print_sum(seq)

    def display_hint(self):
        return 'sum'


class ProdPrinter:

    def __init__(self, val):
        self.seq = val['seq']
        self.overall_coeff = ex_to_number(val['overall_coeff'])

    def should_skip_overall_coeff(self):
        return ClNumber(self.overall_coeff).fixum == 1

    def print_term(self, term):
        base = term['rest']
        exponent = ex_to_number(term['coeff'])
        if ClNumber(exponent).fixnum == 1:
            return '({0})'.format(str(base))
        else:
            str_exponent = ClNumberPrinter(exponent).to_string()
            return '({0})^({1})'.format(str(base), str_exponent)

    def print_seq(self):
        return ' * '.join([self.print_term(term) for term in
                           StdVectorIterator(self.seq)])

    def to_string(self):
        if self.should_skip_overall_coeff:
            return self.print_seq()
        else:
            coeff = ClNumberPrinter(self.overall_coeff).to_string()
            return '{0}*({1})'.format(coeff, self.print_seq())


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
