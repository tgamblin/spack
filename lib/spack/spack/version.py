##############################################################################
# Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This file is part of Spack.
# Created by Todd Gamblin, tgamblin@llnl.gov, All rights reserved.
# LLNL-CODE-647188
#
# For details, see https://github.com/llnl/spack
# Please also see the LICENSE file for our notice and the LGPL.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (as
# published by the Free Software Foundation) version 2.1, February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
# conditions of the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##############################################################################
"""
This module implements Version and version-ish objects.  These are:

Version
  A single version of a package.
VersionRange
  A range of versions of a package.
VersionList
  A list of Versions and VersionRanges.

All of these types support the following operations, which can
be called on any of the types::

  __eq__, __ne__, __lt__, __gt__, __ge__, __le__, __hash__
  __contains__
  satisfies
  could_satisfy
  overlaps
  union
  intersection
  exact
  exactly
  copy

"""
import re
import numbers
from bisect import bisect_left
from functools import wraps

from functools_backport import total_ordering
from spack.util.spack_yaml import syaml_dict

__all__ = ['Version', 'VersionRange', 'VersionList', 'ver']

# Valid version characters
VALID_VERSION = r'[A-Za-z0-9_.-]'


def int_if_int(string):
    """Convert a string to int if possible.  Otherwise, return a string."""
    try:
        return int(string)
    except ValueError:
        return string


def coerce_versions(a, b):
    """
    Convert both a and b to the 'greatest' type between them, in this order:
           Version < VersionRange < VersionList
    This is used to simplify comparison operations below so that we're always
    comparing things that are of the same type.
    """
    order = (Version, VersionRange, VersionList)
    ta, tb = type(a), type(b)

    def check_type(t):
        if t not in order:
            raise TypeError("coerce_versions cannot be called on %s" % t)
    check_type(ta)
    check_type(tb)

    if ta == tb:
        return (a, b)
    elif order.index(ta) > order.index(tb):
        if ta == VersionRange:
            return (a, VersionRange(b, b))
        else:
            return (a, VersionList([b]))
    else:
        if tb == VersionRange:
            return (VersionRange(a, a), b)
        else:
            return (VersionList([a]), b)


def coerced(method):
    """Decorator that ensures that argument types of a method are coerced."""
    @wraps(method)
    def coercing_method(a, b, *args, **kwargs):
        if type(a) == type(b) or a is None or b is None:
            return method(a, b, *args, **kwargs)
        else:
            ca, cb = coerce_versions(a, b)
            return getattr(ca, method.__name__)(cb, *args, **kwargs)
    return coercing_method


def _numeric_lt(self0, other):
    """Compares two versions, knowing they're both numeric"""


def eq_ignore_exact(v1, v2):
    return (v1 is None and v2 is None) or (
        v1 is not None and v2 is not None and v1.version == v2.version)


@total_ordering
class Version(object):
    """Class to represent versions"""

    def __init__(self, string, exact=None):
        # Support copying other versions.
        if isinstance(string, Version):
            self.string = string.string
            self.version = string.version
            self.separators = string.separators
            self._exact = string._exact if exact is None else bool(exact)
            return

        string = str(string)

        if not re.match(VALID_VERSION, string):
            raise ValueError("Bad characters in version string: %s" % string)

        # preserve the original string, but trimmed.
        string = re.sub(r'[^a-zA-Z0-9.]*$', '', string.strip())
        self.string = string

        # Split version into alphabetical and numeric segments
        segment_regex = r'[a-zA-Z]+|[0-9]+'
        segments = re.findall(segment_regex, string)
        self.version = tuple(int_if_int(seg) for seg in segments)

        # Store the separators from the original version string as well.
        # last element of separators is ''
        separators = re.split(segment_regex, string)
        self.separators = tuple(separators[1:-1])

        # Support exact versions literals. Exact versions have a final period.
        self._exact = exact or (separators[-1] == '.')
        if self._exact:
            # strip period from exact version literal
            self.string = self.string.rstrip('.')

    @property
    def dotted(self):
        return '.'.join(str(x) for x in self.version)

    @property
    def underscored(self):
        return '_'.join(str(x) for x in self.version)

    @property
    def dashed(self):
        return '-'.join(str(x) for x in self.version)

    @property
    def joined(self):
        return ''.join(str(x) for x in self.version)

    def up_to(self, index):
        """Return a version string up to the specified component, exclusive.
           e.g., if this is 10.8.2, self.up_to(2) will return '10.8'.
        """
        return '.'.join(str(x) for x in self[:index])

    def lowest(self):
        return self

    def highest(self):
        return self

    def isnumeric(self):
        """Tells if this version is numeric (vs. a non-numeric version).  A
        version will be numeric as long as the first section of it is,
        even if it contains non-numerica portions.

        Some numeric versions:
            1
            1.1
            1.1a
            1.a.1b
        Some non-numeric versions:
            develop
            system
            myfavoritebranch
        """
        return isinstance(self.version[0], numbers.Integral)

    def isdevelop(self):
        """Triggers on the special case of the `@develop` version."""
        return self.string == 'develop'

    @coerced
    def satisfies(self, other):
        """A Version 'satisfies' another if it is at least as specific and has
        a common prefix.  e.g., we want gcc@4.7.3 to satisfy a request for
        gcc@4.7 so that when a user asks to build with gcc@4.7, we can find
        a suitable compiler.

        An *exact* version will ONLY be satisfied with
        """
        if other.exact:
            return self._exact and self.version == other.version
        else:
            nself = len(self.version)
            nother = len(other.version)
            return nother <= nself and self.version[:nother] == other.version

    @coerced
    def could_satisfy(self, other):
        """This version *could* satisfy another if it becomes more specific."""
        return other.satisfies(self)

    def wildcard(self):
        """Create a regex that will match variants of this version string."""
        def a_or_n(seg):
            if type(seg) == int:
                return r'[0-9]+'
            else:
                return r'[a-zA-Z]+'

        version = self.version

        # Use a wildcard for separators, in case a version is written
        # two different ways (e.g., boost writes 1_55_0 and 1.55.0)
        sep_re = '[_.-]'
        separators = ('',) + (sep_re,) * len(self.separators)

        version += (version[-1],) * 2
        separators += (sep_re,) * 2

        segments = [a_or_n(seg) for seg in version]

        wc = segments[0]
        for i in xrange(1, len(separators)):
            wc += '(?:' + separators[i] + segments[i]

        # Add possible alpha or beta indicator at the end of each segemnt
        # We treat these specially b/c they're so common.
        wc += '(?:[a-z]|alpha|beta)?)?' * (len(segments) - 1)
        return wc

    def __iter__(self):
        return iter(self.version)

    def __getitem__(self, idx):
        cls = type(self)
        if isinstance(idx, numbers.Integral):
            return self.version[idx]
        elif isinstance(idx, slice):
            # Currently len(self.separators) == len(self.version) - 1
            extendend_separators = self.separators + ('',)
            string_arg = []
            for token, sep in zip(self.version, extendend_separators)[idx]:
                string_arg.append(str(token))
                string_arg.append(str(sep))
            string_arg.pop()  # We don't need the last separator
            string_arg = ''.join(string_arg)
            return cls(string_arg)
        message = '{cls.__name__} indices must be integers'
        raise TypeError(message.format(cls=cls))

    @property
    def yaml_str(self):
        """In YAML, we keep the final '.' for exact versions."""
        rstring = self.string
        if self._exact:
            rstring += '.'
        return rstring

    def __repr__(self):
        """Representation for the interpreter."""
        return 'Version(%r)' % self.yaml_str

    def __str__(self):
        """String to display to the user."""
        return self.string

    @property
    def exact(self):
        """A Version is exact if it could represent only one version.

        Examples:
          - 2.7 could mean 2.7.12, 2.7.11, 2.7.2.1, etc., so it's not exact.
          - 2.7. means *exactly* 2.7, NOT 2.7.*, so it is exact.

        If this version is exact, this returns it.  Otherwise it returns None.
        """
        return self if self._exact else None

    @property
    def exactly(self):
        """Return a copy of this Version with ``exact`` set to True."""
        if self._exact:
            return self
        else:
            return Version(self, exact=True)

    @property
    def inexactly(self):
        if not self._exact:
            return self
        else:
            return Version(self, exact=False)

    def copy(self):
        return Version(self)

    def _numeric_lt(self, other):
        """Compares two versions, knowing they're both numeric"""
        # Standard comparison of two numeric versions
        for a, b in zip(self.version, other.version):
            if a == b:
                continue
            else:
                # Numbers are always "newer" than letters.
                # This is for consistency with RPM.  See patch
                # #60884 (and details) from bugzilla #50977 in
                # the RPM project at rpm.org.  Or look at
                # rpmvercmp.c if you want to see how this is
                # implemented there.
                if type(a) != type(b):
                    return type(b) == int
                else:
                    return a < b
        # If the common prefix is equal, the one
        # with more segments is bigger.
        return len(self.version) < len(other.version)

    @coerced
    def __lt__(self, other):
        """Version comparison is designed for consistency with the way RPM
           does things.  If you need more complicated versions in installed
           packages, you should override your package's version string to
           express it more sensibly.

           An exact version is treated as less than the identical inexact
           version. i.e.: ``2.7.3.`` < ``2.7.3``
        """
        if other is None:
            return False

        # simple equality test first.
        if self.version == other.version:
            return self._exact and not other._exact

        # TODO: revisit this develop stuff.
        # First priority: anything < develop
        if self.isdevelop():
            return False    # source = develop, it can't be < anything

        # Now we know !sdev
        if other.isdevelop():
            return True    # src < dst

        # now we know neither self nor other isdevelop().

        # Principle: Non-numeric is less than numeric
        # (so numeric will always be preferred by default)
        if self.isnumeric():
            if other.isnumeric():
                return self._numeric_lt(other)
            else:    # self = numeric; other = non-numeric
                # Numeric > Non-numeric (always)
                return False
        else:
            if other.isnumeric():  # self = non-numeric, other = numeric
                # non-numeric < numeric (always)
                return True
            else:  # Both non-numeric
                # Maybe consider other ways to compare here...
                return self.string < other.string

    @coerced
    def __eq__(self, other):
        return (other is not None and
                type(other) == Version and
                self._exact == other._exact and
                self.version == other.version)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.version)

    @coerced
    def __contains__(self, other):
        if other is None:
            return False

        if self._exact:
            return other and other.exact and self.version == other.version
        else:
            return other.version[:len(self.version)] == self.version

    def is_predecessor(self, other):
        """True if the other version is the immediate predecessor of this one.
           That is, NO versions v exist such that:
           (self < v < other and v not in self).
        """
        if len(self.version) != len(other.version):
            return False

        # exact versions always have more specific versions between them
        # and any other version.
        if self.exact:
            return False

        sl = self.version[-1]
        ol = other.version[-1]
        return type(sl) == int and type(ol) == int and (ol - sl == 1)

    def is_successor(self, other):
        return other.is_predecessor(self)

    @coerced
    def overlaps(self, other):
        return self in other or other in self

    @coerced
    def union(self, other):
        if self == other or other in self:
            return self
        elif self in other:
            return other
        else:
            return VersionList([self, other])

    @coerced
    def intersection(self, other):
        if eq_ignore_exact(self, other):
            if self is other or self == other:
                return self
            elif self.exact:
                return self
            else:  # other.exact
                return other
        else:
            return VersionList()


@total_ordering
class VersionRange(object):

    def __init__(self, start, end):
        if isinstance(start, basestring):
            start = Version(start)
        if isinstance(end, basestring):
            end = Version(end)

        self.start = start
        self.end = end
        if start and end and end < start:
            raise ValueError("Invalid Version range: %s" % self)

    def lowest(self):
        return self.start

    def highest(self):
        return self.end

    @coerced
    def __lt__(self, other):
        """Sort VersionRanges lexicographically so that they are ordered first
           by start and then by end.  None denotes an open range, so None in
           the start position is less than everything except None, and None in
           the end position is greater than everything but None.
        """
        if other is None:
            return False

        s, o = self, other
        # as with __eq__, compare tuples here to ignore cases where the
        # lower bound is exact (as it's the same as inexact lower bound)
        if not eq_ignore_exact(s.start, o.start):
            return s.start is None or (
                o.start is not None and s.start < o.start)
        return (s.end != o.end and
                o.end is None or (s.end is not None and s.end < o.end))

    @coerced
    def __eq__(self, other):
        return (other is not None and
                type(other) == VersionRange and
                # Note that (2.7.3.:2.8) == (2.7.3:2.8)
                # b/c 2.7.3. is the smallest possible version in 2.7.3
                eq_ignore_exact(self.start, other.start) and
                # but, (2.7.3:2.8.) is NOT the same as (2.7.3:2.8)
                self.end == other.end)

    def __ne__(self, other):
        return not (self == other)

    @property
    def exact(self):
        """Whether this VersionRange represents an exact version, or could
           potentially represent many versions."""
        if (self.start == self.end) and self.start and self.start.exact:
            return self.start
        else:
            return None

    @property
    def exactly(self):
        """If this VersionRange can represent an exact version, return it."""
        if not eq_ignore_exact(self.start, self.end):
            return None

        if self.start.exact and self.end.exact:
            return self
        else:
            return VersionRange(self.start.exactly, self.end.exactly)

    @property
    def inexactly(self):
        """Return self for inexact ranges. Convert exact ranges to inexact."""
        if (not self.start or not self.start.exact) or (
                not self.end or not self.end.exact):
            return self

        return VersionRange(Version(self.start, exact=False),
                            Version(self.end, exact=False))

    @coerced
    def __contains__(self, other):
        if other is None:
            return False

        in_lower = (self.start == other.start or
                    self.start is None or
                    (other.start is not None and (
                        self.start < other.start or
                        other.start in self.start)))
        if not in_lower:
            return False

        in_upper = (self.end == other.end or
                    self.end is None or
                    (other.end is not None and (
                        self.end > other.end or
                        other.end in self.end)))
        return in_upper

    @coerced
    def could_satisfy(self, other):
        """A VersionRange could satisfy another if some version in this range
        could satisfy some version in other's range.

        So, either:

        a) It overlaps with the other range
        b) The start of this range satisfies the end of the other range.

        This is the same as ``overlaps()``.
        """
        return self.overlaps(other)

    @coerced
    def satisfies(self, other):
        """A VersionRange satisfies another if all versions in the range are in
           the other range.

           This is equivalent to ``self in other``.
        """
        return self in other

    @coerced
    def overlaps(self, other):
        return ((self.start is None or other.end is None or
                 self.start <= other.end or
                 other.end in self.start or self.start in other.end) and
                (other.start is None or self.end is None or
                 other.start <= self.end or
                 other.start in self.end or self.end in other.start))

    @coerced
    def union(self, other):
        if not self.overlaps(other):
            if (self.end is not None and other.start is not None and
                    self.end.is_predecessor(other.start)):
                return VersionRange(self.start, other.end)

            if (other.end is not None and self.start is not None and
                    other.end.is_predecessor(self.start)):
                return VersionRange(other.start, self.end)

            return VersionList([self, other])

        # if we're here, then we know the ranges overlap.
        if self.start is None or other.start is None:
            start = None
        else:
            start = self.start
            # TODO: See note in intersection() about < and in discrepancy.
            if self.start in other.start or other.start < self.start:
                start = other.start

        if self.end is None or other.end is None:
            end = None
        else:
            end = self.end
            # TODO: See note in intersection() about < and in discrepancy.
            if other.end not in self.end:
                if end in other.end or other.end > self.end:
                    end = other.end

        return VersionRange(start, end)

    @coerced
    def intersection(self, other):
        if self.overlaps(other):
            if self.start is None:
                start = other.start
            else:
                start = self.start
                if other.start is not None:
                    if other.start > start or other.start in start:
                        start = other.start

            if self.end is None:
                end = other.end
            else:
                end = self.end
                # Handling end case is tricky:
                #
                #  1. 1.6.5 in 1.6 = True  (1.6.5 is more specific)
                #  2. 1.6 < 1.6.5  = True  (lexicographic)
                #  3. 1.6. < 1.6           (exact less than inexact)
                #
                # Here we test (not end in other.end) to handle #2
                if other.end is not None and end not in other.end:
                    if other.end < end or other.end in end:
                        end = other.end

                # Handle #3: if the ranges only overlap at an exact
                # version, we need to constrain more.
                if eq_ignore_exact(end, start):
                    if start.exact or end.exact:
                        start = end = start.exactly

            # Return a Version if start and end are equal
            if start and end and start == end:
                return start

            # Now we've got a VersionRange.
            return VersionRange(start, end)

        else:
            return VersionList()

    def __hash__(self):
        return hash((self.start, self.end))

    def copy(self):
        return VersionRange(self)

    @property
    def yaml_str(self):
        out = ''
        if self.start:
            out += self.start.yaml_str
        out += ':'
        if self.end:
            out += self.end.yaml_str
        return out

    def __repr__(self):
        return 'VersionRange(%r)' % self.yaml_str

    def __str__(self):
        out = ''
        if self.start:
            out += str(self.start)
        out += ':'
        if self.end:
            out += str(self.end)
        return out


@total_ordering
class VersionList(object):
    """Sorted, non-redundant list of Versions and VersionRanges."""

    def __init__(self, vlist=None):
        self.versions = []
        if vlist is not None:
            if isinstance(vlist, basestring):
                vlist = _string_to_version(vlist)
                if type(vlist) == VersionList:
                    self.versions = vlist.versions
                else:
                    self.versions = [vlist]
            else:
                vlist = list(vlist)
                for v in vlist:
                    self.add(ver(v))

    def add(self, version):
        if type(version) in (Version, VersionRange):
            # This normalizes single-value version ranges.
            if version.exact:
                version = version.exact

            i = bisect_left(self, version)

            while i - 1 >= 0 and version.overlaps(self[i - 1]):
                version = version.union(self[i - 1])
                del self.versions[i - 1]
                i -= 1

            while i < len(self) and version.overlaps(self[i]):
                version = version.union(self[i])
                del self.versions[i]

            self.versions.insert(i, version)

        elif type(version) == VersionList:
            for v in version:
                self.add(v)

        else:
            raise TypeError("Can't add %s to VersionList" % type(version))

    @property
    def concrete(self):
        return self.exact

    @property
    def exact(self):
        """A VersionList is exact if it can represent exactly ONE version."""
        return self[0].exact if len(self) == 1 else None

    @property
    def exactly(self):
        """A VersionList is exact if it can represent exactly ONE version."""
        if len(self) == 1:
            return self if self[0].exact else VersionList([self[0].exactly])
        return None

    @property
    def inexactly(self):
        """Return self for inexact lists. Convert exact lists to inexact."""
        if any(not v.exact for v in self):
            return self
        return VersionList([v.inexactly for v in self])

    def copy(self):
        return VersionList(self)

    def lowest(self):
        """Get the lowest version in the list."""
        if not self:
            return None
        else:
            return self[0].lowest()

    def highest(self):
        """Get the highest version in the list."""
        if not self:
            return None
        else:
            return self[-1].highest()

    @coerced
    def overlaps(self, other):
        if not other or not self:
            return False

        s = o = 0
        while s < len(self) and o < len(other):
            if self[s].overlaps(other[o]):
                return True
            elif self[s] < other[o]:
                s += 1
            else:
                o += 1
        return False

    def to_dict(self):
        """Generate human-readable dict for YAML."""
        if self.exact:
            return syaml_dict([
                ('version', self[0].yaml_str)
            ])
        else:
            return syaml_dict([
                ('versions', self.yaml_str)
            ])

    @staticmethod
    def from_dict(dictionary):
        """Parse dict from to_dict."""
        if 'versions' in dictionary:
            return VersionList(dictionary['versions'])
        elif 'version' in dictionary:
            return VersionList([dictionary['version']])
        else:
            raise ValueError("Dict must have 'version' or 'versions' in it.")

    @coerced
    def satisfies(self, other, strict=True):
        """A VersionList satisfies another if all versions in the list are
           in the other.

           This is equivalent to ``self in other``.
        """
        # TODO: consider factoring out `strict` parameter.
        if not strict:
            return self.could_satisfy(other)

        return self in other

    @coerced
    def could_satisfy(self, other):
        """A VersionList could satisfy another if any version in this
           list is in the other.

           This is equivalent to ``self.overlaps(other)``.
        """
        return self.overlaps(other)

    @coerced
    def update(self, other):
        for v in other.versions:
            self.add(v)

    @coerced
    def union(self, other):
        result = self.copy()
        result.update(other)
        return result

    @coerced
    def intersection(self, other):
        # TODO: make this faster.  This is O(n^2).
        result = VersionList()
        for s in self:
            for o in other:
                result.add(s.intersection(o))
        return result

    @coerced
    def intersect(self, other):
        """Intersect this spec's list with other.

        Return True if the spec changed as a result; False otherwise
        """
        isection = self.intersection(other)
        changed = (isection.versions != self.versions)
        self.versions = isection.versions
        return changed

    @coerced
    def __contains__(self, other):
        if len(self) == 0:
            return False

        for version in other:
            i = bisect_left(self, version)
            if i == 0:
                if version not in self[0]:
                    return False
            elif all(version not in v for v in self[i - 1:]):
                return False

        return True

    def __getitem__(self, index):
        return self.versions[index]

    def __iter__(self):
        return iter(self.versions)

    def __reversed__(self):
        return reversed(self.versions)

    def __len__(self):
        return len(self.versions)

    @coerced
    def __eq__(self, other):
        return other is not None and self.versions == other.versions

    def __ne__(self, other):
        return not (self == other)

    @coerced
    def __lt__(self, other):
        return other is not None and self.versions < other.versions

    def __hash__(self):
        return hash(tuple(self.versions))

    def __str__(self):
        return ",".join(str(v) for v in self.versions)

    @property
    def yaml_str(self):
        return ",".join(v.yaml_str for v in self.versions)

    def __repr__(self):
        return 'VersionList(%r)' % self.yaml_str


def _string_to_version(string):
    """Converts a string to a Version, VersionList, or VersionRange.
       This is private.  Client code should use ver().
    """
    string = string.replace(' ', '')

    if ',' in string:
        return VersionList(string.split(','))

    elif ':' in string:
        s, e = string.split(':')
        start = Version(s) if s else None
        end = Version(e) if e else None
        return VersionRange(start, end)

    else:
        return Version(string)


def ver(obj):
    """Parses a Version, VersionRange, or VersionList from a string
       or list of strings.
    """
    if isinstance(obj, (list, tuple)):
        return VersionList(obj)
    elif isinstance(obj, basestring):
        return _string_to_version(obj)
    elif isinstance(obj, (int, float)):
        return _string_to_version(str(obj))
    elif type(obj) in (Version, VersionRange, VersionList):
        return obj
    else:
        raise TypeError("ver() can't convert %s to version!" % type(obj))
