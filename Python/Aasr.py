#!/usr/bin/python3
"""Azimuth-Ambiguity-to-Signal Ratio
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = ""

import numpy as np
import math

import Constants
import Utils

# -----------------------------------------------------------------------------
# AASR.
def Get_AASR(DopplerBandwidth, Frequencies, NumAmbiguities, FreqSteps, Prf,
    TxAntennaPattern, RxAntennaPattern, TxGain=1, RxGain=1):
    """
    Calculate a number of Azimuth-Ambiguity-to-Signal Ratios with a given
    processed Doppler Bandwidth, radar platform velocity and a frequency steps
    with a given antenna pattern for transmit & rececive , with PRF.
    
    DopplerBandwidth in Hertz
    NumAmbiguities is a scalar
    FreqSteps is a scalar
    Prf in Hertz
    TxAntennaPattern is the magnitude of the TX antenna pattern over the 
        DopplerBandwidth range.
        This is a dictionray of form: dict{position: magnitude}
    RxAntennaPattern is the magnitude of the RX antenna pattern over the 
        DopplerBandwidth range
        This is a dictionray of form: dict{position: magnitude}
    
    returns the AASR as a magnitude value
    """
    _BWDop = DopplerBandwidth
    _NA = NumAmbiguities
    _step = _BWDop / FreqSteps
    
    _denom = 0
    for f in np.arange(start=-_BWDop/2, stop=_BWDop/2, step=_step):
        # find the nearest frequency and magnitude pair for the TX and RX patterns
        _idx, _ = Utils.find_nearest(Frequencies, f)
        _txAzPat = TxAntennaPattern[_idx] * TxGain
        _rxAzPat = RxAntennaPattern[_idx] * RxGain
        _azPat = _txAzPat * _rxAzPat
        _denom += ((math.pow(_azPat, 2)) * _step)

        # _txFreq, _txAzPat = Utils.find_nearest(list(TxAntennaPattern.keys()), f)
        # _rxFreq, _rxAzPat = Utils.find_nearest(list(RxAntennaPattern.keys()), f)

        # # _txAzPat = Get_AntennaPatternWithDopplerFrequency(f, txDimension, velocity, txGain)
        # # _rxAzPat = Get_AntennaPatternWithDopplerFrequency(f, rxDimension, velocity, rxGain)
        # # _denom += ((math.pow(_azPat, 2)) * _step)
        # _denom += ((_txAzPat * _rxAzPat) * _step)
        
    _num = 0
    _azPat = 0
    for m in range(-_NA, _NA+1, 1):
        if m != 0:
            for f in np.arange(start=-_BWDop/2, stop=_BWDop/2, step=_step):
                # find the nearest frequency and magnitude pair for the TX and RX patterns
                _idx, _ = Utils.find_nearest(Frequencies, f + (m*Prf))
                _txAzPat = TxAntennaPattern[_idx] * TxGain
                _rxAzPat = RxAntennaPattern[_idx] * RxGain
                _azPat = _txAzPat * _rxAzPat
                _num += ((math.pow(_azPat, 2)) * _step)

                # _txFreq, _txAzPat = Utils.find_nearest(list(TxAntennaPattern.keys()), f + (m*Prf))
                # _rxFreq, _rxAzPat = Utils.find_nearest(list(RxAntennaPattern.keys()), f + (m*Prf))

                # # _txAzPat = Get_AntennaPatternWithDopplerFrequency((f + (m*prf)), txDimension, velocity, txGain)
                # # _rxAzPat = Get_AntennaPatternWithDopplerFrequency((f + (m*prf)), rxDimension, velocity, rxGain)
                # # #_num += ((math.pow(_azPat, 2)) * _step)
                # _num += ((_txAzPat * _rxAzPat) * _step)
    
    return (_num / _denom)

