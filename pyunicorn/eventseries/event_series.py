#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of pyunicorn.
# Copyright (C) 2008--2019 Jonathan F. Donges and pyunicorn authors
# URL: <http://www.pik-potsdam.de/members/donges/software>
# License: BSD (3-clause)
#
# Please acknowledge and cite the use of this software and its authors
# when results are used in publications or published elsewhere.
#
# You can use the following reference:
# J.F. Donges, J. Heitzig, B. Beronov, M. Wiedermann, J. Runge, Q.-Y. Feng,
# L. Tupikina, V. Stolbova, R.V. Donner, N. Marwan, H.A. Dijkstra,
# and J. Kurths, "Unified functional network and nonlinear time series analysis
# for complex systems science: The pyunicorn package"

"""
Provides class for event series analysis, namely event synchronization and event
coincidence analysis.

"""

#
# Imports
#

import numpy as np

from .. import cached_const

import warnings


class EventSeries:


    def __init__(self, eventmatrix, taumax):
        """
        Initialize an instance of EventSeries.

        Format of eventmatrix:
        An eventmatrix is a 2D numpy array with the first dimension covering
        the timesteps and the second dimensions covering the variables. Each
        variable at a specific timestep is either '1' if an event occured or
        '0' if it did not, e.g. for 3 variables with 10 timesteps the
        eventmatrix could look like

            array([[0, 1, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [1, 0, 1],
                   [0, 1, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [0, 1, 0],
                   [0, 0, 0]])

        Important note!!!
        Due to normalization issues the event synchronization matrices should
        only be used if the number of events (i.e the number 1s) in each
        variable is identical.

        :type eventmatrix: 2D Numpy array [time, variables]
        :arg eventmatrix: Event series array
        :arg int taumax: Maximum dynamical delay
        """

        self.__T = eventmatrix.shape[0]
        self.__N = eventmatrix.shape[1]
        self.__eventmatrix = eventmatrix
        self.taumax = taumax

        # Dictionary for chached constants
        self.cache = {'base': {}}
        """(dict) cache of re-usable computation results"""

        # Check for right input format
        if len(np.unique(eventmatrix)) != 2 or not (
                np.unique(eventmatrix) == np.array([0, 1])).all():
            raise ValueError("Eventmatrix not in correct format")

        # Print warning if number of events is not identical for all variables
        NrOfEvs = np.sum(eventmatrix, axis=0)
        if not (NrOfEvs == NrOfEvs[0]).all():
            warnings.warn("Data does not contain equal number of events")


    @staticmethod
    def _eventsync(EventSeriesX, EventSeriesY, taumax=None, taulag=0):
        """
        Calculates the directed event synchronization from two event series X
        and Y.

        :type EventSeriesX: 1D Numpy array
        :arg EventSeriesX: Event series containing '0's and '1's
        :type EventSeriesY: 1D Numpy array
        :arg EventSeriesY: Event series containing '0's and '1's
        :rtype: list
        :return: [Event synchronization XY, Event synchronization YX]
        """

        if taumax is None:
            taumax = np.inf

        # Get time indices (type boolean or simple '0's and '1's)
        ex = np.array(np.where(EventSeriesX), dtype=np.int8)
        ey = np.array(np.where(EventSeriesY), dtype=np.int8)
        # Number of events
        lx = ex.shape[1]
        ly = ey.shape[1]
        if lx == 0 or ly == 0:              # Division by zero in output
            return np.nan, np.nan
        if lx in [1, 2] or ly in [1, 2]:    # Too few events to calculate
            return 0., 0.
        # Array of distances
        dstxy2 = 2 * (np.repeat(ex[:, 1:-1].T, ly-2, axis=1)
                      - np.repeat(ey[:, 1:-1], lx-2, axis=0))
        # Dynamical delay
        diffx = np.diff(ex)
        diffy = np.diff(ey)
        diffxmin = np.minimum(diffx[:, 1:], diffx[:, :-1])
        diffymin = np.minimum(diffy[:, 1:], diffy[:, :-1])
        tau2 = np.minimum(np.repeat(diffxmin.T, ly-2, axis=1),
                          np.repeat(diffymin, lx-2, axis=0))
        tau2 = np.minimum(tau2, 2 * taumax)
        # Count equal time events and synchronised events
        eqtime = dstxy2.size - np.count_nonzero(dstxy2)

        # Calculate boolean matrices of coincidences
        Axy = (dstxy2 > 0) * (dstxy2 <= tau2)
        Ayx = (dstxy2 < 0) * (dstxy2 >= -tau2)

        # Loop over coincidences and determine number of double counts
        # by checking at least one event of the pair is also coincided
        # in other direction
        countxydouble = countyxdouble = 0

        for i, j in np.transpose(np.where(Axy)):
            countxydouble += np.any(Ayx[i, :]) or np.any(Ayx[:, j])
        for i, j in np.transpose(np.where(Ayx)):
            countyxdouble += np.any(Axy[i, :]) or np.any(Axy[:, j])

        # Calculate counting quantities and subtract half of double countings
        countxy = np.sum(Axy) + 0.5 * eqtime - 0.5 * countxydouble
        countyx = np.sum(Ayx) + 0.5 * eqtime - 0.5 * countyxdouble

        norm = np.sqrt((lx-2) * (ly-2))
        return countxy / norm, countyx / norm


    @staticmethod
    def _eca(EventSeriesX, EventSeriesY, delT, tau=0, ts1=None, ts2=None):
        """
        Event coincidence analysis:
        Returns the precursor and trigger coincidence rates of two event series
        X and Y.

        :type EventSeriesX: 1D Numpy array
        :arg EventSeriesX: Event series containing '0's and '1's
        :type EventSeriesY: 1D Numpy array
        :arg EventSeriesY: Event series containing '0's and '1's
        :arg delT: coincidence interval width
        :arg int tau: lag parameter
        :rtype: list
        :return: [Precursor coincidence rate XY, Trigger coincidence rate XY,
              Precursor coincidence rate YX, Trigger coincidence rate YX]
        """

        # Count events that cannot be coincided due to tau and delT
        if not (tau == 0 and delT == 0):
            # Start of EventSeriesX
            n11 = np.count_nonzero(EventSeriesX[:tau + delT])
            # End of EventSeriesX
            n12 = np.count_nonzero(EventSeriesX[-(tau + delT):])
            # Start of EventSeriesY
            n21 = np.count_nonzero(EventSeriesY[:tau + delT])
            # End of EventSeriesY
            n22 = np.count_nonzero(EventSeriesY[-(tau + delT):])
        else:
            # Instantaneous coincidence
            n11, n12, n21, n22 = 0, 0, 0, 0
        # Get time indices
        if ts1 is None:
            e1 = np.where(EventSeriesX)[0]
        else:
            e1 = ts1[EventSeriesX]
        if ts2 is None:
            e2 = np.where(EventSeriesY)[0]
        else:
            e2 = ts2[EventSeriesY]
        del EventSeriesX, EventSeriesY, ts1, ts2
        # Number of events
        l1 = len(e1)
        l2 = len(e2)
        # Array of all interevent distances
        dst = (np.array([e1] * l2).T - np.array([e2] * l1))
        # Count coincidences with array slicing
        prec12 = np.count_nonzero(np.any(((dst - tau >= 0)
                                          * (dst - tau <= delT))[n11:, :],
                                         axis=1))
        trig12 = np.count_nonzero(np.any(((dst - tau >= 0)
                                          * (dst - tau <= delT))
                                         [:, :dst.shape[1] - n22],
                                         axis=0))
        prec21 = np.count_nonzero(np.any(((-dst - tau >= 0)
                                          * (-dst - tau <= delT))[:, n21:],
                                         axis=0))
        trig21 = np.count_nonzero(np.any(((-dst - tau >= 0)
                                          * (-dst - tau <= delT))
                                         [:dst.shape[0] - n12, :],
                                         axis=1))
        # Normalisation and output
        return (np.float32(prec12) / (l1 - n11), np.float32(trig12) / (l2 - n22),
                np.float32(prec21) / (l2 - n21), np.float32(trig21) / (l1 - n12))
