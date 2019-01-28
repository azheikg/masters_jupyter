#!/usr/bin/python3
"""Range-Ambiguity-to-Signal Ratio
"""

__author__ = "Kevin Gema"
__copyright__ = "Copyright 2018, Kevin Gema"
__credits__ = ["Kevin Gema"]
__license__ = "GPL"
__version__ = "1.0.2"
__maintainer__ = ""

import numpy as np
import math

import Constants
import Utils
import SphericalEarthModel as SEM

# -----------------------------------------------------------------------------
# Ambiguous slant range
def AmbiguousSlantRange(desiredSlantRange, multiple, prf):
    """
    Calculates the ambiguous slant range for a given desired slant range, 
    multiple of PRF (ambiguity index) and prf.
    
    desiredSlantRange in meters
    multiple is a scalar
    prf in Hertz
    """
    return (desiredSlantRange + (multiple * (Constants.SPEED_OF_LIGHT / (2*prf))))


# -----------------------------------------------------------------------------
# RASR.
def Get_RasrMain(SlantRange, IncidenceAngle,
    TxAntennaPattern, RxAntennaPattern):
    num = np.power(SlantRange, 3) * np.sin(IncidenceAngle)
    txAzPat = TxAntennaPattern
    rxAzPat = RxAntennaPattern
    pat = np.abs(txAzPat * rxAzPat)
    # pat = txAzPat * rxAzPat
    antPat = np.power(pat, 2)
    # antPat = pat
    return num / antPat

def CalRasrMainAntennaPatterns_alt(rxAnt, txAnt, wavelength, numScanAngles, numTxAngles, minLookAngle, 
    maxLookAngle, numSamples, txGain=-1, rxGain=-1):
    """
    Calculate the antenna field patterns over the TX/RX beamwidths
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

            # get the receciver field pattern for this scan angle, over the calculated angle range
            if numScanAngles > 1:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + j*txAngleInc
            else:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + (rxElBW/2)  + j*txAngleInc
            rxPat = rxAnt.elevation_field_pattern(wavelength, anglesRad, rxSteer, 1, True) * rxGain
            RxAntPatterns.append(rxPat)

            # get the transmitter field pattern for this scan angle, over the calculated angle range
            txSteer = j*txAngleInc
            txPat = txAnt.elevation_field_pattern(wavelength, anglesRad, txSteer, 1, True) * txGain
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

def CalRasrMainAntennaPatterns_rect(rxAnt, txAnt, wavelength, numScanAngles, numTxAngles, minLookAngle, 
    maxLookAngle, numSamples, txGain=-1, rxGain=-1):
    """
    Calculate the antenna field patterns over the TX/RX beamwidths
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

            # get the receciver field pattern for this scan angle, over the calculated angle range
            if numScanAngles > 1:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + j*txAngleInc
            else:
                rxSteer = (-(txElBW/2) + i*rxAngleInc) + (rxElBW/2)  + j*txAngleInc
            rxPat = rxAnt.elevation_field_pattern(wavelength, anglesRad, rxSteer) * rxGain
            RxAntPatterns.append(rxPat)

            # get the transmitter field pattern for this scan angle, over the calculated angle range
            txSteer = j*txAngleInc
            txPat = txAnt.elevation_field_pattern(wavelength, anglesRad, txSteer) * txGain
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

def CalcAmbiguousAntennaPatterns(rxAnt, txAnt, slantRange, numNearAmb, numFarAmb, altitude, 
    antennaOffset, wavelength, prf, txGain=-1, rxGain=-1):
    """Calculate the slant ranges, incidence angles and antenna patterns for
    the ambiuous signals from the number of near and far ambiguities"""
    SrAmb_Arr = []
    IaAmb_Arr = []
    TxAmb_Arr = []
    RxAmb_Arr = []

    for j in range(-numNearAmb, numFarAmb+1, 1):
        if j != 0:
            # for each ambiguity, calculate the slant range the incidence angle and the look angle
            # then calculate the antenna pattern angle from this look angle, then calculate the 
            # antenna patterns
            # save all these values for the summation calculation
            sr_amb_arr = slantRange + (Constants.SPEED_OF_LIGHT/2) * (j/prf)
            
            # # where the ambiguous slant range is smaller than the altitude,
            # # replace with the altitude, just to compensate for overflows which
            # # impacts the visualisation of the RASR graph
            # sr_amb_arr[sr_amb_arr < altitude] = 0

            ia_amb_arr = SEM.IncidenceAngleFromSlantRange(altitude, sr_amb_arr)
            la_amb_arr = SEM.LookAngle(altitude, ia_amb_arr)
            ant_amb_angle_arr = la_amb_arr - Utils.deg2rad(antennaOffset)

            if txGain == -1:
                # get the antenna gain from the antenna model
                txPat_amb = txAnt.elevation_power_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * txAnt.gain(wavelength)
            else:
                txPat_amb = txAnt.elevation_power_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * txGain

            if rxGain == -1:
                # get the antenna gain from the antenna model
                rxPat_amb = rxAnt.elevation_power_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * rxAnt.gain(wavelength)
            else:
                rxPat_amb = rxAnt.elevation_power_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * rxGain

            SrAmb_Arr.append(sr_amb_arr)
            IaAmb_Arr.append(ia_amb_arr)
            TxAmb_Arr.append(txPat_amb)
            RxAmb_Arr.append(rxPat_amb)
    # convert lists to arrays
    SrAmb_Arr = np.asarray(SrAmb_Arr)
    IaAmb_Arr = np.asarray(IaAmb_Arr)
    TxAmb_Arr = np.asarray(TxAmb_Arr)
    RxAmb_Arr = np.asarray(RxAmb_Arr)

    return TxAmb_Arr, RxAmb_Arr, SrAmb_Arr, IaAmb_Arr

def CalcAmbiguousAntennaPatterns_alt(rxAnt, txAnt, numScanAngles, slantRange, numNearAmb, numFarAmb, altitude, 
    antennaOffset, wavelength, prf, txGain=-1, rxGain=-1):
    """Calculate the slant ranges, incidence angles and antenna patterns for
    the ambiuous signals from the number of near and far ambiguities"""

    SR_Amb_Scan_Arr = []
    IA_Amb_Scan_Arr = []
    TX_Amb_Scan_Arr = []
    RX_Amb_Scan_Arr = []

    for i in range(numScanAngles):
        SrAmb_Arr = []
        IaAmb_Arr = []
        TxAmb_Arr = []
        RxAmb_Arr = []

        for j in range(-numNearAmb, numFarAmb+1, 1):
            if j != 0:
                # for each ambiguity, calculate the slant range the incidence angle and the look angle
                # then calculate the antenna pattern angle from this look angle, then calculate the 
                # antenna patterns
                # save all these values for the summation calculation
                sr_amb_arr = slantRange[i] + (Constants.SPEED_OF_LIGHT/2) * (j/prf)
                
                # # where the ambiguous slant range is smaller than the altitude,
                # # replace with the altitude, just to compensate for overflows which
                # # impacts the visualisation of the RASR graph
                # sr_amb_arr[sr_amb_arr < altitude] = 0

                ia_amb_arr = SEM.IncidenceAngleFromSlantRange(altitude, sr_amb_arr)
                la_amb_arr = SEM.LookAngle(altitude, ia_amb_arr)
                ant_amb_angle_arr = la_amb_arr - Utils.deg2rad(antennaOffset)

                if txGain == -1:
                    # get the antenna gain from the antenna model
                    txPat_amb = txAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * txAnt.gain(wavelength)
                else:
                    txPat_amb = txAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * txGain

                if rxGain == -1:
                    # get the antenna gain from the antenna model
                    rxPat_amb = rxAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * rxAnt.gain(wavelength)
                else:
                    rxPat_amb = rxAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0, 1, True) * rxGain

                SrAmb_Arr.append(sr_amb_arr)
                IaAmb_Arr.append(ia_amb_arr)
                TxAmb_Arr.append(txPat_amb)
                RxAmb_Arr.append(rxPat_amb)

        # append to the  Scan arrays
        SR_Amb_Scan_Arr.append(SrAmb_Arr)
        IA_Amb_Scan_Arr.append(IaAmb_Arr)
        TX_Amb_Scan_Arr.append(TxAmb_Arr)
        RX_Amb_Scan_Arr.append(RxAmb_Arr)

    # convert lists to arrays
    SR_Amb_Scan_Arr = np.asarray(SR_Amb_Scan_Arr)
    IA_Amb_Scan_Arr = np.asarray(IA_Amb_Scan_Arr)
    TX_Amb_Scan_Arr = np.asarray(TX_Amb_Scan_Arr)
    RX_Amb_Scan_Arr = np.asarray(RX_Amb_Scan_Arr)

    return TX_Amb_Scan_Arr, RX_Amb_Scan_Arr, SR_Amb_Scan_Arr, IA_Amb_Scan_Arr

def CalcAmbiguousAntennaPatterns_rect(rxAnt, txAnt, numScanAngles, slantRange, numNearAmb, numFarAmb, altitude, 
    antennaOffset, wavelength, prf, txGain=-1, rxGain=-1):
    """Calculate the slant ranges, incidence angles and antenna patterns for
    the ambiuous signals from the number of near and far ambiguities"""

    SR_Amb_Scan_Arr = []
    IA_Amb_Scan_Arr = []
    TX_Amb_Scan_Arr = []
    RX_Amb_Scan_Arr = []

    for i in range(numScanAngles):
        SrAmb_Arr = []
        IaAmb_Arr = []
        TxAmb_Arr = []
        RxAmb_Arr = []

        for j in range(-numNearAmb, numFarAmb+1, 1):
            if j != 0:
                # for each ambiguity, calculate the slant range the incidence angle and the look angle
                # then calculate the antenna pattern angle from this look angle, then calculate the 
                # antenna patterns
                # save all these values for the summation calculation
                sr_amb_arr = slantRange[i] + (Constants.SPEED_OF_LIGHT/2) * (j/prf)
                
                # # where the ambiguous slant range is smaller than the altitude,
                # # replace with the altitude, just to compensate for overflows which
                # # impacts the visualisation of the RASR graph
                # sr_amb_arr[sr_amb_arr < altitude] = 0

                ia_amb_arr = SEM.IncidenceAngleFromSlantRange(altitude, sr_amb_arr)
                la_amb_arr = SEM.LookAngle(altitude, ia_amb_arr)
                ant_amb_angle_arr = la_amb_arr - Utils.deg2rad(antennaOffset)

                if txGain == -1:
                    # get the antenna gain from the antenna model
                    txPat_amb = txAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0) * txAnt.gain(wavelength)
                else:
                    txPat_amb = txAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0) * txGain

                if rxGain == -1:
                    # get the antenna gain from the antenna model
                    rxPat_amb = rxAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0) * rxAnt.gain(wavelength)
                else:
                    rxPat_amb = rxAnt.elevation_field_pattern(wavelength, ant_amb_angle_arr, 0) * rxGain

                SrAmb_Arr.append(sr_amb_arr)
                IaAmb_Arr.append(ia_amb_arr)
                TxAmb_Arr.append(txPat_amb)
                RxAmb_Arr.append(rxPat_amb)

        # append to the  Scan arrays
        SR_Amb_Scan_Arr.append(SrAmb_Arr)
        IA_Amb_Scan_Arr.append(IaAmb_Arr)
        TX_Amb_Scan_Arr.append(TxAmb_Arr)
        RX_Amb_Scan_Arr.append(RxAmb_Arr)

    # convert lists to arrays
    SR_Amb_Scan_Arr = np.asarray(SR_Amb_Scan_Arr)
    IA_Amb_Scan_Arr = np.asarray(IA_Amb_Scan_Arr)
    TX_Amb_Scan_Arr = np.asarray(TX_Amb_Scan_Arr)
    RX_Amb_Scan_Arr = np.asarray(RX_Amb_Scan_Arr)

    return TX_Amb_Scan_Arr, RX_Amb_Scan_Arr, SR_Amb_Scan_Arr, IA_Amb_Scan_Arr


def RasrAmb(TxPat, RxPat, SR, IA):
    """
    Calculate the ambiguity inner sum
    """
    antPat = np.power(np.abs(TxPat * RxPat), 2)
    return antPat / (np.power(SR, 3) * np.sin(IA))