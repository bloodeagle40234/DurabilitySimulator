#!/usr/bin/python
# Copyright (c) 2014 Kota Tsuyuzaki

"""
An Implements of Durability Calculator for OpenStack Swift.

To build a very large scale object storage system, we have to consider
about the suitable design such as processing HW, network architecture,
# of devices and system configuration.

For these design, sometimes we need to know how durability is calculated
to confirm that our storage system is really durable or not. However it could
be difficalut since the components of the storage system would be complicated
for beginers and there is a lot of reference papers for the durability
calculation.

This script helps them to calculate how their own storage system is durable.
Its calculation is built based on a paper[1]. There is a similar idea called
"redundancy set" with the swift's partition. We can assumes a redundancy set
as a partition for durability calculation. Please see the paper in detail.

[1]: "Reliability Mechanisms for Very Large Storage Systems"
      http://www.ssrc.ucsc.edu/Papers/xin-mss03.pdf
"""

import numpy as np
import math
from decimal import Decimal
from optparse import OptionParser


def generate_malkov_inv_matrix(k, m, u, v, print_title=False):
    """
    This function returns inverse of M (i.e. I-Q) matrix
    instead of Q for ensuring precision. M is a limitation of
    state transition matrix among temprary states.

    :param: k: # of data fragments.
    :param: m: # of parity fragments.
    """

    # k must be an integer more than 0
    if k < 1:
        raise ValueError

    # n: # of all fragments
    n = k + m

    # At Reaplica Model (case k==1)
    # we have to all
    loosable = m + 1 if k == 1 else m

    tmp = []
    for x in xrange(loosable):
        row = []
        row_title = []
        for y in xrange(loosable):
            # skip unreachable statement
            # TODO: append 0 on this state
            if abs(x - y) > 1:
                row_title.append('0')
                row.append(0)
                continue
            if x == y:
                ux = n - x
                cell = Decimal(1) - (Decimal(ux) * Decimal(u) +
                                     Decimal(x) * Decimal(v))
                cell = float(Decimal(1) - Decimal(cell))
                row.append(cell)
                row_title.append('1-(%du + %dv)' % (ux, x))
            elif x > y:
                row.append(x * v * -1)
                row_title.append('%dv' % x)
            elif y > x:
                val = n - x
                row.append((val) * u * -1)
                row_title.append('%du' % val)
        if print_title:
            print row_title
        tmp.append(row)
    return np.matrix(tmp).getI()


def generate_state_matrix(k, m):
    """
    Generates m rows, 1 column matrix filled by 1.
    """
    loosable = m + 1 if k == 1 else m
    state_matrix = np.matrix([[1] for x in xrange(loosable)])
    return state_matrix


if __name__ == '__main__':
    """
    Replica Model (2 replicas):
        k = 1, m = 1
    Replica Model (3 replicas):
        k = 1, m = 2
    EC Model (6 data fragments and 3 parity fragments):
        k = 6, m = 3

    Note that EC model could allow to lose m fragments but Replica Model
    could allow m+1 fragments. To simplify, this script specializes Replica
    Model case (i.e. k == 1) to calculate the duability.
    """

    parser = OptionParser()

    parser.add_option(
        '-k', '--data-num', dest='data_num', default=1,
        type=int, help='# of data fragments. Default is 1. (3 Replica Model)')
    parser.add_option(
        '-m', '--parity-num', dest='parity_num', default=2, type=int,
        help='# of parity fragments. Default is 2. (3 Replica Model)')
    parser.add_option(
        '-p', '--parition_num', dest='part_power', default=18,
        type=int, help='# of partition power. Default is 18.')
    parser.add_option(
        '-u', '--afr', dest='afr', default=8.3/10**6, type=float,
        help='Disk AFR from catalog (or actual) spec (per hour)')
    # TODO: add repair rate calculation description with network design
    parser.add_option(
        '-v', '--repair-rate', dest='repair_rate', default=1.66,
        type=float, help='Repair Rate')
    parser.add_option(
        '--vorbose', dest='vorbose', default=False,
        type=float, help='show detail calculation method')

    # parse args and set the values
    (options, args) = parser.parse_args()
    data_num = options.data_num
    parity_num = options.parity_num
    rs_num = 2 ** options.part_power
    u = options.afr
    # u = 1.4 / 10 ** 8
    v = options.repair_rate
    # v = 0.0277
    vorbose = options.vorbose

    print '########## Parameters for Calculation ##########'
    print '# of data fragments: %d' % data_num
    print '# of parity fragments: %d' % parity_num
    print '# of partition in Swift: %d' % rs_num
    print 'AFR: %f' % u
    print 'RepairRate: %f' % v
    print '############# Calculation Results ##############'
    # generate matrix from malcov model
    if vorbose:
        print 'Malkov matrix'
    matrix = generate_malkov_inv_matrix(data_num, parity_num, u, v, vorbose)
    state_matrix = generate_state_matrix(data_num, parity_num)

    mttdl_rs = matrix * state_matrix
    mttdl = mttdl_rs / rs_num
    durability = 1 - (1 / mttdl.getA1()[0])
    try:
        nines = int(math.log(1 - durability, 0.1))
    except ValueError:
        nines = 'Over flow (more than 15)'

    if vorbose:
        print 'mttdl_rs: \n', mttdl_rs
        print 'mttdl: \n', mttdl
    print 'Storage Efficiency (NOTE: no replica -> 1.0): %s' % \
        (float(data_num + parity_num) / float(data_num))
    print 'Swift Durability : %s' % durability
    print '# of Nines (LOG): %s' % nines
    print '################################################'
