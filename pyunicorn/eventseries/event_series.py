#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of pyunicorn.
# Copyright (C) 2008--2020 Jonathan F. Donges and pyunicorn authors
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
coincidence analysis. The data format event matrix has to be passed in order to
instantiate an object of the class which may contain an arbitray number of event
series. However, the methods event synchronization and event coincidence analysis
may be called without instantiating an object of the class.

"""

#
# Imports
#

import numpy as np

from .. import cached_const

import warnings


class EventSeries:

    def __init__(self, eventmatrix):
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
        """

        self.__T = eventmatrix.shape[0]
        self.__N = eventmatrix.shape[1]
        self.__eventmatrix = eventmatrix

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

    def __str__(self):
        """
        Return a string representation of the EventSeries object.
        """
        return ('EventSynchronization: %i variables, %i timesteps'
                % (self.__N, self.__T))

    @staticmethod
    def event_sync(eventseriesx, eventseriesy, ts1=None, ts2=None, taumax=np.inf, step=0, windowsize=0):
        """
        Calculates the directed event synchronization from two event series X
        and Y.

        :type eventseriesx: 1D Numpy array
        :arg eventseriesx: Event series containing '0's and '1's
        :type eventseriesy: 1D Numpy array
        :arg eventseriesy: Event series containing '0's and '1's
        :type ts1: 1D Numpy array
        :arg ts1: Event time array containing time points when events of event series 1 occur, not obligatory
        :type ts2: 1D Numpy array
        :arg ts2: Event time array containing time points when events of event series 2 occur, not obligatory
        :type taumax:
        :arg taumax: maximum distance of two events to be counted as synchronous
        :type step:
        :arg step:
        :type windowsize:
        :arg windowsize:
        :rtype: list
        :return: [Event synchronization XY, Event synchronization YX]
        """
        # Get time indices
        if ts1 is None:
            e1 = np.where(eventseriesx)[0]
        else:
            e1 = ts1[eventseriesx]
        if ts2 is None:
            e2 = np.where(eventseriesy)[0]
        else:
            e2 = ts2[eventseriesy]
        # del es1, es2, ts1 , ts2
        # Number of events
        l1 = len(e1)
        l2 = len(e2)
        if l1 == 0 or l2 == 0:  # Division by zero in output
            return np.nan, np.nan
        if l1 in [1, 2] or l2 in [1, 2]:  # Too few events to calculate
            return 0., 0.


        dst12 = (np.array([e1[1:-1]] * (l2 - 2), dtype='int32').T - np.array([e2[1:-1]] * (l1 - 2), dtype='int32'))
        # Dynamical delay
        diff1 = np.diff(e1)
        diff2 = np.diff(e2)
        diff1min = np.minimum(diff1[1:], diff1[:-1])
        diff2min = np.minimum(diff2[1:], diff2[:-1])
        tau = 0.5 * np.minimum(np.array([diff1min] * (l2 - 2), dtype='float32').T,
                               np.array([diff2min] * (l1 - 2), dtype='float32'))
        efftau = np.minimum(tau, taumax)
        # del diff1 , diff2 , diff1min, diff2min, tau
        # Count equal time events and synchronised events
        eqtime = 0.5 * np.sum(dst12 == 0)
        count12 = np.sum((dst12 > 0) * (dst12 <= efftau)) + eqtime
        count21 = np.sum((dst12 < 0) * (-dst12 <= efftau)) + eqtime
        norm = np.sqrt((l1 - 2) * (l2 - 2))
        # del dst12
        return np.float32(count12 / norm), np.float32(count21 / norm)

    @cached_const('base', 'directedES')
    def directed_event_sync(self, ):
        """
        Returns the NxN matrix of the directed event synchronization measure.
        The entry [i, j] denotes the directed event synchronization from
        variable j to variable i.
        """
        eventmatrix = self.__eventmatrix
        res = np.ones((self.__N, self.__N)) * np.inf

        for i in range(0, self.__N):
            for j in range(i+1, self.__N):
                res[i, j], res[j, i] = self.eventSync(eventmatrix[:, i],
                                                       eventmatrix[:, j], )
        return res

    def symmetric_event_sync(self):
        """
        Returns the NxN matrix of the undirected or symmetrix event
        synchronization measure. It is obtained by the sum of both directed
        versions.
        """
        directed = self.directed_event_sync()
        return directed + directed.T

    def antisymmetric_event_sync(self):
        """
        Returns the NxN matrix of the antisymmetric synchronization measure.
        It is obtained by the difference of both directed versions.
        """
        directed = self.directedES()
        return directed - directed.T

    @staticmethod
    def event_coincidence_analysis(eventseriesx, eventseriesy, ts1=None, ts2=None, tau=0, delT=3, step=0, windowsize=0):
        """
         Event coincidence analysis:
         Returns the precursor and trigger coincidence rates of two event series
         X and Y.

         :type eventseriesx: 1D Numpy array
         :arg eventseriesx: Event series containing '0's and '1's
         :type eventseriesy: 1D Numpy array
         :arg eventseriesy: Event series containing '0's and '1's
         :arg delT: coincidence interval width
         :arg int tau: lag parameter
         :type ts1: 1D Numpy array
         :arg ts1: Event time array containing time points when events of event series 1 occur, not obligatory
         :type ts2: 1D Numpy array
         :arg ts2: Event time array containing time points when events of event series 2 occur, not obligatory
         :type step:
         :arg step:
         :type windowsize:
         :arg windowsize:
         :rtype list
         :return [Precursor coincidence rate XY, Trigger coincidence rate XY,
               Precursor coincidence rate YX, Trigger coincidence rate YX]
         """
        # Count events that cannot be coincided due to tau and delT
        if not (tau == 0 and delT == 0):
            n11 = np.count_nonzero(eventseriesx[:tau + delT])  # Start of es1
            n12 = np.count_nonzero(eventseriesx[-(tau + delT):])  # End of es1
            n21 = np.count_nonzero(eventseriesy[:tau + delT])  # Start of es2
            n22 = np.count_nonzero(eventseriesy[-(tau + delT):])  # End of es2
        else:
            n11, n12, n21, n22 = 0, 0, 0, 0  # Instantaneous coincidence
        # Get time indices
        if ts1 is None:
            e1 = np.where(eventseriesx)[0]
        else:
            e1 = ts1[eventseriesx]
        if ts2 is None:
            e2 = np.where(eventseriesy)[0]
        else:
            e2 = ts2[eventseriesy]
        # del es1, es2, ts1, ts2

        # Number of events
        l1 = len(e1)
        l2 = len(e2)
        # Array of all interevent distances
        dst = (np.array([e1] * l2).T - np.array([e2] * l1))
        # print "shape of interevent distances matrix:", dst.shape
        # Count coincidences with array slicing
        # print dst - tau > 0# , (dst - tau <= delT)
        prec12 = np.count_nonzero(
            np.any(((dst - tau >= 0) * (dst - tau <= delT))[n11:, :], axis=1))  # CHECK [n11:,:] = [n11:,] = [n11:] ??
        trig12 = np.count_nonzero(np.any(((dst - tau >= 0) *
                                          (dst - tau <= delT))[:, :dst.shape[1] - n22], axis=0))
        prec21 = np.count_nonzero(np.any(((-dst - tau >= 0) *
                                          (-dst - tau <= delT))[:, n21:], axis=0))
        trig21 = np.count_nonzero(np.any(((-dst - tau >= 0) *
                                          (-dst - tau <= delT))[:dst.shape[0] - n12, :], axis=1))
        # Normalisation and output
        del dst
        return (np.float32(prec12) / (l1 - n11), np.float32(trig12) / (l2 - n22), np.float32(prec21) / (l2 - n21),
                np.float32(trig21) / (l1 - n12))
