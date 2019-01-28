#!/usr/bin/python3
"""Noise-Equivalent Sigma Zero
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.2"
__maintainer__ = ""

import numpy as np
import math

import SphericalEarthModel as SEM
import Constants
import Utils

# -----------------------------------------------------------------------------
# NESZ
def NESZ(PlatformVelocity, SquintAngle, SlantRange, NPSD, SysLoss, TxPowerAvg, 
    TxAntennaGain, RxAntennaGain, wavelength, SlantRangeResolution, 
    GrazingAngle):
    """
    Calculate the NESZ for a given slant range, radar platform velocity, 
    squint angle Noise Power Spectral Density (NPSD), System Losses (SysLoss),
    Average transmitter power, transmit and receive antenna gains, 
    wavelength, slant range resolution and grazing angle (GA)
    
    PlatformVelocity in m/s
    SquintAngle in radians
    SlantRange in meters
    NPSD is a scalar
    SysLoss in dB
    TxPowerAvg in Watts
    TxAntennaGain is a scalar
    RxAntennaGain is a scalar
    wavelength in meters
    SlantRangeResolution in meters
    GrazingAngle in radians
    """
    R = SlantRange
    V = PlatformVelocity
    GA = GrazingAngle
    PAvg = TxPowerAvg
    Gt = TxAntennaGain
    Gr = RxAntennaGain
    SrRes = SlantRangeResolution

    num = 2 * (np.power((4 * np.pi * R), 3)) * NPSD * math.pow(10, SysLoss/10) \
            * np.sin(SquintAngle) * V * np.cos(GA)

    denom = PAvg * Gt * Gr * math.pow(wavelength, 3) * SrRes

    return (num / denom)

def NESZ2(Pn, Prf, Ia, Naz, wavelength, Pavg, Gt, Gr, AzRes, R, V, CrxEl, CtxEl, CrxAz, CtxAz):
    # _naz = wavelength * R * Prf / (2*AzRes * V)
    # print(_naz)
    _num1 = 2*np.power((4*np.pi),3) * Pn * Prf * Naz * np.sin(Ia)
    _denom1 = Constants.SPEED_OF_LIGHT * (wavelength**2) * Pavg * Gt*Gr * AzRes

    _num2 = 1
    _C2wayEl = CrxEl * CtxEl

    _innerSum = 0
    # print(Naz)
    for i in range(np.trunc(Naz)):
        _C2wayAz = CrxAz * CtxAz
        _inner = _C2wayAz / np.powerr(R, 2)
        _innerSum += _inner
    
    if _innerSum == 0:
        _innerSum = 1

    _denom2 = np.power(np.abs(_C2wayEl * _innerSum), 2)

    return ((_num1/_denom1) * (_num2/_denom2))


def NESZ_BRET(R, V, Ia, Lsa, T, Brg, F, L, Ptx, Gtx, Grx, wavelength, Pl, Prf):
    """
    Using the formula for NESZ presented in 'Bistatic Radar: Emerging Technology', p117
    """
    _F = 10**(F/10)
    _L = 10**(L/10)
    _num = (4**4)*np.power((np.pi*R), 3)* V*np.sin(Ia - Lsa)*Constants.BOLTZMANS_CONSTANT*T*Brg*_F*_L
    _denom = Ptx*Gtx*Grx*(wavelength**3)*Constants.SPEED_OF_LIGHT*Pl*Prf
    return _num/_denom


# -----------------------------------------------------------------------------
def CalNeszAntennaPatterns(rxAnt, txAnt, numScanAngles, txBeamSteer=False, minLookAngle=None, maxLookAngle=None):
    """
    Calculate the antenna patterns over the TX/RX beamwidths
    """
    RxAntPatterns = []
    TxAntPatterns = []
    rxtxAnglesRad = []
    rxPeaks = {}

    txElBW = txAnt.elevation_beamwidth(Wavelength)
    rxElBW = rxAnt.elevation_beamwidth(Wavelength)
    txGain = txAnt.gain(Wavelength)
    rxGain = rxAnt.gain(Wavelength)
    
    # calculate the increment between different stripmap sessions.
    # this is effectively the overlap between swaths
    # the -1 is to ensure there is a overlap at the 2 beam peaks
    if txBeamSteer:
        assert minLookAngle is not None, "Must specify min/max look angles when doing beam steering"
        assert maxLookAngle is not None, "Must specify min/max look angles when doing beam steering"
        #rxtxAngleInc = ((LookAngleRangeRad.max() - LookAngleRangeRad.min()) / (numScanAngles -1)) 
        rxtxAngleInc = ((maxLookAngle - minLookAngle) / (numScanAngles -1)) 
    else:
        rxtxAngleInc = (txElBW) / (numScanAngles -1) # the -1 is to ensure there is a overlap at the 2 beam peaks

    # calculate the antenna patterns for each of the different scan angles
    # each pattern is calculated over the range 'anglesRad'
    for i in range(numScanAngles):
        # get the angle range for this scan angle
        anglesRad = np.linspace((-(txElBW/2) + i*rxtxAngleInc) - (rxElBW/2), 
                                (-(txElBW/2) + i*rxtxAngleInc) + (rxElBW/2), 
                                NumSamples)
        rxtxAnglesRad.append(anglesRad)

        # get the receiver power pattern for this scan angle, over the calculated angle range
        rxPat = rxAnt.elevation_power_pattern(Wavelength, anglesRad, (-(txElBW/2) + i*rxtxAngleInc), 1, True) * rxGain
        RxAntPatterns.append(rxPat)

        # get the transmitter power pattern for this scan angle, over the calculated angle range
        if txBeamSteer:
            txPat = txAnt.elevation_power_pattern(Wavelength, anglesRad, (-(txElBW/2) + i*rxtxAngleInc), 1, True) * txGain
        else:
            txPat = txAnt.elevation_power_pattern(Wavelength, anglesRad, 0, 1, True) * txGain
        TxAntPatterns.append(txPat)

        # keep track of the peaks of the scanned RX patterns
        for idx, agl in enumerate(anglesRad):
            # to angles in radians are very small
            # so keep track of a larger version, and compensate later
            angle = int(agl*1000)  
            if angle in rxPeaks:
                if rxPeaks[angle] < rxPat[idx]:
                    rxPeaks[angle] = rxPat[idx]
            else:
                rxPeaks[angle] = rxPat[idx]

    # Get the maximum RX pattern values and the angles where they occur
    RxAntPatMaxAngles = np.fromiter(rxPeaks.keys(), dtype=float)
    RxAntPatMaxValues = np.fromiter(rxPeaks.values(), dtype=float)
    
    return RxAntPatterns,TxAntPatterns, rxtxAnglesRad, rxPeaks, RxAntPatMaxAngles, RxAntPatMaxValues

def CalNeszAntennaPatterns_alt(rxAnt, txAnt, wavelength, numScanAngles, numTxAngles, minLookAngle, 
    maxLookAngle, numSamples, txGain=-1, rxGain=-1):
    """
    Calculate the antenna patterns over the TX/RX beamwidths
    """
    RxAntPatterns = []
    TxAntPatterns = []
    rxtxAnglesRad = []
    rxPeaks = {}

    txElBW = txAnt.elevation_beamwidth(wavelength)
    rxElBW = rxAnt.elevation_beamwidth(wavelength)
    if txGain == -1:
        # get the antenna gain from the antenna model
        txGain = txAnt.gain(wavelength)
    if rxGain == -1:
        # get the antenna gain from the antenna model
        rxGain = rxAnt.gain(wavelength)
    
    # divide the look angle range between the number of TX angles
    txAngleInc = ((maxLookAngle - minLookAngle) / (numTxAngles)) 

    # for each of the TX angles, there is 'numScanAngles' RX angles
    if numScanAngles > 1:
        rxAngleInc = (txElBW) / (numScanAngles-1) # the -1 is to ensure there is a overlap at the 2 beam peaks
    else:
        rxAngleInc = txAngleInc + (rxElBW/2)

    # for each of the TX angles, scan the RX beam over these angles
    for j in range(numTxAngles):
        for i in range(numScanAngles):
            # calc the rx angles for each TX beam
            if numScanAngles > 1:
                anglesRad = np.linspace((-(txElBW/2) + i*rxAngleInc) - (rxElBW/2) + j*txAngleInc,
                                          (-(txElBW/2) + i*rxAngleInc) + (rxElBW/2) + j*txAngleInc,
                                          numSamples)
            else:
                anglesRad = np.linspace((-(txElBW/2) + i*rxAngleInc) + j*txAngleInc,
                                          (-(txElBW/2) + i*rxAngleInc) + (rxElBW) + j*txAngleInc,
                                          numSamples)
            rxtxAnglesRad.append(anglesRad)

            # get the receciver power pattern for this scan angle, over the calculated angle range
            if numScanAngles > 1:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + j*txAngleInc
            else:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + (rxElBW/2)  + j*txAngleInc
            rxPat = rxAnt.elevation_power_pattern(wavelength, anglesRad, rxSteer, 1, True) * rxGain
            RxAntPatterns.append(rxPat)

            # get the transmitter power pattern for this scan angle, over the calculated angle range
            txSteer = j*txAngleInc
            txPat = txAnt.elevation_power_pattern(wavelength, anglesRad, txSteer, 1, True) * txGain
            TxAntPatterns.append(txPat)

            # keep track of the peaks of the scanned RX patterns
            for idx, agl in enumerate(anglesRad):
                # to angles in radians are very small
                # so keep track of a larger version, and compensate later
                angle = int(agl*1000)  
                if angle in rxPeaks:
                    if rxPeaks[angle] < rxPat[idx]:
                        rxPeaks[angle] = rxPat[idx]
                else:
                    rxPeaks[angle] = rxPat[idx]

            # Get the maximum RX pattern values and the angles where they occur
            RxAntPatMaxAngles = np.fromiter(rxPeaks.keys(), dtype=float)
            RxAntPatMaxValues = np.fromiter(rxPeaks.values(), dtype=float)
    
    return RxAntPatterns,TxAntPatterns, rxtxAnglesRad, rxPeaks, RxAntPatMaxAngles, RxAntPatMaxValues

def CalNeszAntennaPatterns_rect(rxAnt, txAnt, wavelength, numScanAngles, numTxAngles, minLookAngle, 
    maxLookAngle, numSamples, txGain=-1, rxGain=-1):
    """
    Calculate the antenna patterns over the TX/RX beamwidths for rectangular apertures
    """
    RxAntPatterns = []
    TxAntPatterns = []
    rxtxAnglesRad = []
    rxPeaks = {}

    txElBW = txAnt.elevation_beamwidth(wavelength)
    rxElBW = rxAnt.elevation_beamwidth(wavelength)
    if txGain == -1:
        # get the antenna gain from the antenna model
        txGain = txAnt.gain(wavelength)
    if rxGain == -1:
        # get the antenna gain from the antenna model
        rxGain = rxAnt.gain(wavelength)
    
    # divide the look angle range between the number of TX angles
    txAngleInc = ((maxLookAngle - minLookAngle) / (numTxAngles)) 

    # for each of the TX angles, there is 'numScanAngles' RX angles
    if numScanAngles > 1:
        rxAngleInc = (txElBW) / (numScanAngles-1) # the -1 is to ensure there is a overlap at the 2 beam peaks
    else:
        rxAngleInc = txAngleInc + (rxElBW/2)

    # for each of the TX angles, scan the RX beam over these angles
    for j in range(numTxAngles):
        for i in range(numScanAngles):
            # calc the rx angles for each TX beam
            if numScanAngles > 1:
                anglesRad = np.linspace((-(txElBW/2) + i*rxAngleInc) - (rxElBW/2) + j*txAngleInc,
                                          (-(txElBW/2) + i*rxAngleInc) + (rxElBW/2) + j*txAngleInc,
                                          numSamples)
            else:
                anglesRad = np.linspace((-(txElBW/2) + i*rxAngleInc) + j*txAngleInc,
                                          (-(txElBW/2) + i*rxAngleInc) + (rxElBW) + j*txAngleInc,
                                          numSamples)
            rxtxAnglesRad.append(anglesRad)

            # get the receciver power pattern for this scan angle, over the calculated angle range
            if numScanAngles > 1:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + j*txAngleInc
            else:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + (rxElBW/2)  + j*txAngleInc
            rxPat = rxAnt.elevation_power_pattern(wavelength, anglesRad, rxSteer) * rxGain
            RxAntPatterns.append(rxPat)

            # get the transmitter power pattern for this scan angle, over the calculated angle range
            txSteer = j*txAngleInc
            txPat = txAnt.elevation_power_pattern(wavelength, anglesRad, txSteer) * txGain
            TxAntPatterns.append(txPat)

            # keep track of the peaks of the scanned RX patterns
            for idx, agl in enumerate(anglesRad):
                # to angles in radians are very small
                # so keep track of a larger version, and compensate later
                angle = int(agl*1000)  
                if angle in rxPeaks:
                    if rxPeaks[angle] < rxPat[idx]:
                        rxPeaks[angle] = rxPat[idx]
                else:
                    rxPeaks[angle] = rxPat[idx]

            # Get the maximum RX pattern values and the angles where they occur
            RxAntPatMaxAngles = np.fromiter(rxPeaks.keys(), dtype=float)
            RxAntPatMaxValues = np.fromiter(rxPeaks.values(), dtype=float)
    
    return RxAntPatterns,TxAntPatterns, rxtxAnglesRad, rxPeaks, RxAntPatMaxAngles, RxAntPatMaxValues


def CalcAnglesAndSlantRanges(numRxScanAngles, altitude, angles, minLookAngle):
    """
    Calculate the look angles, incidence angles and slant ranges
    over the antenna pattern beamwidth range
    """
    IA_Arr = []
    SR_Arr = []
    for i in range(numRxScanAngles):
        ia_arr = SEM.IncidenceAngle(altitude, angles[i] + minLookAngle)
        sr_arr = SEM.SlantRange(altitude, ia_arr, angles[i] + minLookAngle)
        IA_Arr.append(ia_arr)
        SR_Arr.append(sr_arr)
    return IA_Arr, SR_Arr

def CalcGroundRange(numRxScanAngles, altitude, slantRanges):
    """
    calculate the core angles subtended by the slant ranges, and from these the ground ranges
    """
    GR_Arr = [] 
    for i in range(numRxScanAngles):
        core_angle_arr = SEM.CoreAngle(altitude, slantRanges[i])
        gr_arr = SEM.GroundRangeFromCoreAngle(core_angle_arr)
        GR_Arr.append(gr_arr)
    return GR_Arr

def CalcNesz(numRxScanAngles, altitude, chirpBandwidth, ias, srs, txPat, rxPat, wl, txPeakPow, pl, prf):
    """
    Calulate the value of NESZ of the look angles
    """
    neszPeaks = {}
    NESZ_Arr = []
    for i in range(numRxScanAngles):
        nesz_arr = NESZ_BRET(R = srs[i],
                            V = SEM.PlatformVelocity(altitude),
                            Ia = ias[i],
                            Lsa = 0,
                            T = Constants.STANDARD_TEMPERATURE,
                            Brg = chirpBandwidth,
                            F = 4.3,
                            L = 4.1,
                            Ptx = txPeakPow,
                            Gtx = txPat[i],
                            Grx = rxPat[i],
                            wavelength = wl,
                            Pl = pl,
                            Prf = prf)
        NESZ_Arr.append(nesz_arr)

        # keep track of the NESZ peak values
        for idx, ia in enumerate(ias[i]):
            incidenceAngle = int(ia*1000)
            if incidenceAngle in neszPeaks:
                if neszPeaks[incidenceAngle] > nesz_arr[idx]:
                    neszPeaks[incidenceAngle] = nesz_arr[idx]
            else:
                neszPeaks[incidenceAngle] = nesz_arr[idx]

    # convert to numpy array
    NESZ_Arr = np.asarray(NESZ_Arr)
    # get the power values
    Nesz_dB_Arr = 10*np.log10(NESZ_Arr)

    # Get the maximum NESZ values and the angles where they occur
    Nesz_max_incidence_angles = np.fromiter(neszPeaks.keys(), dtype=float)
    Nesz_max_values = np.fromiter(neszPeaks.values(), dtype=float)
    
    return NESZ_Arr, Nesz_max_incidence_angles, Nesz_max_values
