package rs2d.sequence.common;

import rs2d.commons.log.Log;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;


/**
 * The Current Design does not allow one to use TofSat and SatBand
 * at same time, But one may improve it afterwards (if necessary)
 *
 */


public class TofSat extends SatBand {
    private SeqPrep parent;

    protected boolean isAttAuto;

    protected enum UP implements GeneratorParamEnum {
        TOF2D_ENABLED,
        TOF2D_SB_THICKNESS,
        TOF2D_SB_TX_SHAPE,
        TOF2D_SB_GRAMP_SP,
        TOF2D_SB_DISTANCE_FROM_SLICE,
        TOF2D_SB_POSITION;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
    }

    public TofSat(SeqPrep parent) {
        super(parent);
        this.parent = parent;
    }

    @Override
    public void init() throws Exception {
        super.init();
        isTofBandEnabled = parent.getBoolean(CommonUP.MULTI_PLANAR_EXCITATION) && parent.getBoolean(UP.TOF2D_ENABLED); // TOF not allowed in 3D
        parent.getParam(UP.TOF2D_ENABLED).setValue(isTofBandEnabled);

        if (isTofBandEnabled) { //TofBandEnabled has a higher priority than SatBandEnabled
            isSatBandEnabled = false;
            if (parent.hasParam(SatBand.UP.SATBAND_ENABLED))
                parent.getParam(SatBand.UP.SATBAND_ENABLED).setValue(isSatBandEnabled);
        }
        nb_planar_excitation = (parent.getBoolean(CommonUP.MULTI_PLANAR_EXCITATION) ? parent.getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_3D) : 1);
    }

    @Override
    public void prep() throws Exception {
        prepPulse();
        prepGrad();
        prepPulseComp(isTofBandEnabled);
    }

    @Override
    public boolean isEnabled() {
        return isTofBandEnabled;
    }

    protected void prepGradTable() {
        if (isTofBandEnabled) {
            double satband_distance_from_slice = parent.getDouble(UP.TOF2D_SB_DISTANCE_FROM_SLICE);

            double off_center_slice_pos = satband_distance_from_slice + satband_tof_thickness / 2.0; // sat band cranial from voxel
            double off_center_slice_neg = -off_center_slice_pos;  // caudal
            double off_center_slice = 0;
            if ("BELOW THE SLICE".equalsIgnoreCase(parent.getText(UP.TOF2D_SB_POSITION))) {
                off_center_slice = off_center_slice_neg;
            } else if ("ABOVE THE SLICE".equalsIgnoreCase(parent.getText(UP.TOF2D_SB_POSITION))) {
                off_center_slice = off_center_slice_pos;
            }
            double frequency_offset_sat_slice = -grad_amp_satband_mTpm * off_center_slice * (GradientMath.GAMMA / parent.nucleus.getRatio());

            System.out.println("(parent.getBoolean(CommonUP.MULTI_PLANAR_EXCITATION) ? parent.getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_3D) : 1) "+(parent.getBoolean(CommonUP.MULTI_PLANAR_EXCITATION) ? parent.getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_3D) : 1));
            for (int k = 0; k < nb_planar_excitation; k++) {
                offsetFreqSBTable[k] = (parent.pulseTX.getFrequencyOffset(k) * grad_amp_satband_mTpm / parent.gradSlice.getAmplitude_mTpm()) + frequency_offset_sat_slice;
//                    System.out.println("frequency_offset_tof2d[k]  " + offsetFreqSBTable[k]);
            }

            gradAmpSBSliceTable[0] = grad_amp_satband;
            gradAmpSBPhaseTable[0] = 0;
            gradAmpSBReadTable[0] = 0;
            gradAmpSBSliceSpoilerTable[0] = 0;
            gradAmpSBPhaseSpoilerTable[0] = grad_amp_sat_spoiler;
            gradAmpSBReadSpoilerTable[0] = grad_amp_sat_spoiler;
        }
    }

}

