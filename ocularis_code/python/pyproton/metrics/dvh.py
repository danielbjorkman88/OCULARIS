
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import numpy as np
from numba import jit


@jit(nopython=True)
def _build_dvh(dose_tensor: np.ndarray,
               d_max: float,
               mask: np.ndarray,
               delta_d: float):
    # Make sure there is a bin also for the maximum dose value.
    bins = np.arange(0., d_max+delta_d, delta_d)
    cnt = np.zeros_like(bins)

    for i in range(dose_tensor.shape[0]):
        for j in range(dose_tensor.shape[1]):
            for k in range(dose_tensor.shape[2]):
                # Only add the dose to the histogram if it is part of the
                # structure of interest.
                if mask[i, j, k]:
                    ind = int(dose_tensor[i, j, k] // delta_d)
                    cnt[ind] += 1

    # So far we just have a histogram. Now we want to make it cumulative,
    # i. e. the bin for the lowest dose should also contain the counts of
    # the higher dose bins.
    cnt = np.cumsum(cnt[::-1])[::-1]

    return bins, cnt


class DVH:

    def __init__(self, dose: np.ndarray, mask: np.ndarray, delta_d: float = 0.005):
        '''
        Spacing in Gy.
        '''

        assert dose.shape == mask.shape

        self._dose_tensor = dose
        self._mask = mask
        self._delta_d = delta_d

        # These two arrays contain the DVH curve at discrete points.
        # Between the points, the curve is interpolated linearly.
        n_voxels = np.sum(mask)

        bins, cnt = _build_dvh(self._dose_tensor,
                               np.max(self._dose_tensor),
                               self._mask,
                               self._delta_d)
        self._doses = bins
        self._volume_fractions = cnt / n_voxels

    def V(self, dose):
        '''
        Returns the volume fraction receiving at least dose_fraction [Gy] dose.

        Linear interpolation is used.

        Ex.: V(2.4) returns V2.4Gy dose.

        Parameters
        ----------
        dose: array-like
            A float or a numpy array or list of evaluation positions.
            Dose in Gy at which to evaluate the DVH.
        '''
        return np.interp(dose, self._doses, self._volume_fractions)

    def D(self, volume_fraction):
        '''
        Returns the dose of the hottest volume_fraction voxels in Gy.

        Linear interpolation is used.

        Ex.: D(0.02) returns D2 dose.

        Parameters
        ----------
        volume_fraction: array-like
            A float or a numpy array or list of evaluation positions.
        '''
        return np.interp(volume_fraction,
                         self._volume_fractions[::-1],
                         self._doses[::-1])
