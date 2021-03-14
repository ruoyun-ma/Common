package rs2d.sequence.common;

import rs2d.commons.log.Log;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.table.Order;

import java.util.Arrays;
import java.util.List;


public class SatBand implements ModelInterface {
    private SeqPrep parent;
    protected boolean isSatBandEnabled;
    protected boolean isTofBandEnabled = false;

    protected RFPulse pulseTXSatBand;
    protected boolean isAttAuto;
    private boolean isSBorTSEnabled;
    protected int[] position_sli_ph_rea = new int[6];
    protected int nb_satband;
    protected int nb_planar_excitation; // we cannot use SeqPrep.nb_planar_excitation, since it has not been initialized yet

    double[] offsetFreqSBTable;
    double[] gradAmpSBSliceTable;
    double[] gradAmpSBPhaseTable;
    double[] gradAmpSBReadTable;
    double[] gradAmpSBSliceSpoilerTable;
    double[] gradAmpSBPhaseSpoilerTable;
    double[] gradAmpSBReadSpoilerTable;

    protected enum UP implements GeneratorParamEnum {
        SATBAND_ENABLED,
        SATBAND_TX_SHAPE,
        SATBAND_TX_AMP,
        SATBAND_T1,
        SATBAND_TAU,
        SATBAND_THICKNESS,
        SATBAND_ORIENTATION,
        SATBAND_DISTANCE_FROM_FOV,
        SATBAND_GRAD_AMP_SPOILER;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Enable_sb,
        Tx_amp_sb,
        Tx_phase_sb,
        Time_tx_sb,
        Tx_shape_sb,
        Tx_shape_phase_sb,
        Freq_offset_tx_sb,
        Freq_offset_tx_sb_prep,
        Freq_offset_tx_sb_comp,

        Time_grad_ramp_sb,
        Time_grad_sb,

        Grad_amp_sb_phase,
        Grad_amp_sb_read,
        Grad_amp_sb_slice,
        Grad_amp_sb_phase_spoiler,
        Grad_amp_sb_read_spoiler,
        Grad_amp_sb_slice_spoiler,
        Grad_shape_rise_up,
        Grad_shape_rise_down;
    }

    public SatBand(SeqPrep parent) {
        this.parent = parent;
    }

    @Override
    public void init() throws Exception {
        isAttAuto = parent.getBoolean(CommonUP.TX_AMP_ATT_AUTO);
        isSatBandEnabled = parent.getBoolean(UP.SATBAND_ENABLED);
        position_sli_ph_rea = satBandPrep();
    }

    @Override
    public void initFinal() throws Exception {
        isSBorTSEnabled = isSatBandEnabled || isTofBandEnabled;
        parent.set(SP.Enable_sb, isSBorTSEnabled);
        nb_satband = 1;

        if (isSatBandEnabled) {
            nb_satband = (int) Arrays.stream(position_sli_ph_rea).filter(item -> item == 1).count();
        } else if (parent.models.get("TofSat") != null && parent.models.get("TofSat").isEnabled()) {
            nb_satband = nb_planar_excitation;
        }
        // bug on the loop //TODO: XG: JR, Could you please check it? what is the bug on the loop?
        //nb_satband = 1;

        pulseTXSatBand = RFPulse.createRFPulse(parent.getSequence(), CommonSP.Tx_att, SP.Tx_amp_sb, SP.Tx_phase_sb, SP.Time_tx_sb, SP.Tx_shape_sb, SP.Tx_shape_phase_sb, SP.Freq_offset_tx_sb);
        pulseTXSatBand.setShape(parent.getText(UP.SATBAND_TX_SHAPE), parent.nb_shape_points, "90 degree");

        if (isTofBandEnabled) { //TODO: better merge TOF2D_SB_TX_SHAPE with SATBAND_TX_SHAPE one day
            pulseTXSatBand.setShape(parent.getText(TofSat.UP.TOF2D_SB_TX_SHAPE), parent.nb_shape_points, "90 degree");
        }
    }

    @Override
    public void prep() throws Exception {
        prepPulse();
        prepGrad();
        prepPulseComp(isSatBandEnabled);
    }

    @Override
    public void prepFinal() {
        if (isAttAuto) {
            parent.getParam(UP.SATBAND_TX_AMP).setValue(pulseTXSatBand.getAmp180()); //TODO: XG: JR, please check here, some sequence use pulseTXSatBand.getAmp90() instead
        } else
            pulseTXSatBand.setAmp(parent.getDouble(UP.SATBAND_TX_AMP));
    }

    @Override
    public RFPulse getRfPulses() {
        return pulseTXSatBand;
    }

    @Override
    public boolean isEnabled() {
        return isSatBandEnabled;
    }

    protected void prepPulse() {
        // SAT BAND RF pulse
        setSeqParamTime();

        if (isAttAuto)
            clcPulse();
    }

    protected void prepGrad() {
        offsetFreqSBTable = new double[nb_satband];
        gradAmpSBSliceTable = new double[nb_satband];
        gradAmpSBPhaseTable = new double[nb_satband];
        gradAmpSBReadTable = new double[nb_satband];
        gradAmpSBSliceSpoilerTable = new double[nb_satband];
        gradAmpSBPhaseSpoilerTable = new double[nb_satband];
        gradAmpSBReadSpoilerTable = new double[nb_satband];

        prepGradTable();

        // Apply values ot Gradient
        Gradient gradSatBandSlice = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_sb_slice, SP.Time_tx_sb, SP.Grad_shape_rise_up, SP.Grad_shape_rise_down, SP.Time_grad_ramp_sb, parent.nucleus);
        gradSatBandSlice.setAmplitude(gradAmpSBSliceTable);
        gradSatBandSlice.applyAmplitude(isSatBandEnabled ? Order.LoopB : Order.FourLoop);

        Gradient gradSatBandPhase = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_sb_phase, SP.Time_tx_sb, SP.Grad_shape_rise_up, SP.Grad_shape_rise_down, SP.Time_grad_ramp_sb, parent.nucleus);
        gradSatBandPhase.setAmplitude(gradAmpSBPhaseTable);
        gradSatBandPhase.applyAmplitude(isSatBandEnabled ? Order.LoopB : Order.FourLoop);

        Gradient gradSatBandRead = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_sb_read, SP.Time_tx_sb, SP.Grad_shape_rise_up, SP.Grad_shape_rise_down, SP.Time_grad_ramp_sb, parent.nucleus);
        gradSatBandRead.setAmplitude(gradAmpSBReadTable);
        gradSatBandRead.applyAmplitude(isSatBandEnabled ? Order.LoopB : Order.FourLoop);

        Gradient gradSatBandSpoilerSlice = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_sb_slice_spoiler, SP.Time_grad_sb, SP.Grad_shape_rise_up, SP.Grad_shape_rise_down, SP.Time_grad_ramp_sb, parent.nucleus);
        gradSatBandSpoilerSlice.setAmplitude(gradAmpSBSliceSpoilerTable);
        gradSatBandSpoilerSlice.applyAmplitude(isSatBandEnabled ? Order.LoopB : Order.FourLoop);

        Gradient gradSatBandSpoilerPhase = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_sb_phase_spoiler, SP.Time_grad_sb, SP.Grad_shape_rise_up, SP.Grad_shape_rise_down, SP.Time_grad_ramp_sb, parent.nucleus);
        gradSatBandSpoilerPhase.setAmplitude(gradAmpSBPhaseSpoilerTable);
        gradSatBandSpoilerPhase.applyAmplitude(isSatBandEnabled ? Order.LoopB : Order.FourLoop);

        Gradient gradSatBandSpoilerRead = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_sb_read_spoiler, SP.Time_grad_sb, SP.Grad_shape_rise_up, SP.Grad_shape_rise_down, SP.Time_grad_ramp_sb, parent.nucleus);
        gradSatBandSpoilerRead.setAmplitude(gradAmpSBReadSpoilerTable);
        gradSatBandSpoilerRead.applyAmplitude(isSatBandEnabled ? Order.LoopB : Order.FourLoop);

    }

    protected void prepGradTable() {
        double tx_length_sb = 0.0;
        if (parent.hasParam(CommonUP.TX_LENGTH_90)) {
            tx_length_sb = isSBorTSEnabled ? parent.getDouble(CommonUP.TX_LENGTH_90) : parent.minInstructionDelay;
        } else if (parent.hasParam(CommonUP.TX_LENGTH)) {
            tx_length_sb = isSBorTSEnabled ? parent.getDouble(CommonUP.TX_LENGTH) : parent.minInstructionDelay;
        } else {
            Log.error(getClass(), "User Param TX_LENGTH_90 or TX_LENGTH does not exist");
        }

        double tx_bandwidth_factor_sb = 0.0;
        double grad_amp_sat_spoiler = 0.0;
        double satband_thickness = 0.0;
        if (isSatBandEnabled) {
            tx_bandwidth_factor_sb = parent.getTx_bandwidth_factor(UP.SATBAND_TX_SHAPE, CommonUP.TX_BANDWIDTH_FACTOR, CommonUP.TX_BANDWIDTH_FACTOR_3D);
            grad_amp_sat_spoiler = parent.getDouble(UP.SATBAND_GRAD_AMP_SPOILER);
            satband_thickness = parent.getDouble(UP.SATBAND_THICKNESS);
        } else if (isTofBandEnabled) {
            tx_bandwidth_factor_sb = parent.getTx_bandwidth_factor(TofSat.UP.TOF2D_SB_TX_SHAPE, CommonUP.TX_BANDWIDTH_FACTOR, CommonUP.TX_BANDWIDTH_FACTOR_3D);
            grad_amp_sat_spoiler = parent.getDouble(TofSat.UP.TOF2D_SB_GRAMP_SP);
            satband_thickness = parent.getDouble(TofSat.UP.TOF2D_SB_THICKNESS);
        } else {
            Log.error(getClass(), "User Param SATBAND_TX_SHAPE or SATBAND_GRAD_AMP_SPOILER or SATBAND_THICKNESS and TOF2D_SB_TX_SHAPE " +
                    "or TOF2D_SB_GRAMP_SP or TOF2D_SB_THICKNESS do not exist");
        }

        double tx_bandwidth_sb = tx_bandwidth_factor_sb / tx_length_sb;
        double grad_amp_satband = 0;
        double grad_amp_satband_mTpm = 0;

        if (isSBorTSEnabled) {
            Gradient gradSB = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_sb_read, SP.Time_tx_sb, CommonSP.Grad_shape_rise_up, CommonSP.Grad_shape_rise_down, SP.Time_grad_ramp_sb, parent.nucleus);
            if (!gradSB.prepareSliceSelection(tx_bandwidth_sb, satband_thickness)) {
                double satband_thickness_mod = gradSB.getSliceThickness();
                parent.notifyOutOfRangeParam(UP.SATBAND_THICKNESS, satband_thickness_mod, ((NumberParam) parent.getParam(UP.SATBAND_THICKNESS)).getMaxValue(), "Pulse length too short to reach this satband slice thickness");
                satband_thickness = satband_thickness_mod;
            }
            grad_amp_satband = gradSB.getAmplitude();
            grad_amp_satband_mTpm = gradSB.getAmplitude_mTpm();
        }

        // Calculate allocate Gradient amplitude, spoiler, and freq offset
        double satband_distance_from_fov = parent.getDouble(UP.SATBAND_DISTANCE_FROM_FOV);

        if (isSatBandEnabled) {
            //position_sli_ph_rea = satBandPrep();
            int n = 0;
            if (position_sli_ph_rea[0] == 1) { // SB  in slice position Sup
                gradAmpSBSliceTable[n] = grad_amp_satband;
                gradAmpSBPhaseTable[n] = 0;
                gradAmpSBReadTable[n] = 0;
                gradAmpSBSliceSpoilerTable[n] = 0;
                gradAmpSBPhaseSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBReadSpoilerTable[n] = grad_amp_sat_spoiler;
                double off_center_pos = parent.getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_3D) + parent.getDouble(CommonUP.FIELD_OF_VIEW_3D) / 2.0 + satband_distance_from_fov + satband_thickness / 2.0;
                offsetFreqSBTable[n] = new RFPulse().calculateOffsetFreq(grad_amp_satband_mTpm, off_center_pos);
                n += 1;
            }
            if (position_sli_ph_rea[1] == 1) {
                gradAmpSBSliceTable[n] = grad_amp_satband;
                gradAmpSBPhaseTable[n] = 0;
                gradAmpSBReadTable[n] = 0;
                gradAmpSBSliceSpoilerTable[n] = 0;
                gradAmpSBPhaseSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBReadSpoilerTable[n] = grad_amp_sat_spoiler;
                double off_center_neg = parent.getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_3D) - parent.getDouble(CommonUP.FIELD_OF_VIEW_3D) / 2.0 + satband_distance_from_fov + satband_thickness / 2.0;
                offsetFreqSBTable[n] = new RFPulse().calculateOffsetFreq(grad_amp_satband_mTpm, off_center_neg);
                n += 1;
            }
            if (position_sli_ph_rea[2] == 1) {
                gradAmpSBSliceTable[n] = 0;
                gradAmpSBPhaseTable[n] = grad_amp_satband;
                gradAmpSBReadTable[n] = 0;
                gradAmpSBSliceSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBPhaseSpoilerTable[n] = 0;
                gradAmpSBReadSpoilerTable[n] = grad_amp_sat_spoiler;
                double off_center_pos = parent.getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_2D) + parent.getDouble(CommonUP.FIELD_OF_VIEW_PHASE) / 2.0 + satband_distance_from_fov + satband_thickness / 2.0;
                offsetFreqSBTable[n] = new RFPulse().calculateOffsetFreq(grad_amp_satband_mTpm, off_center_pos);
                n += 1;
            }
            if (position_sli_ph_rea[3] == 1) {
                gradAmpSBSliceTable[n] = 0;
                gradAmpSBPhaseTable[n] = grad_amp_satband;
                gradAmpSBReadTable[n] = 0;
                gradAmpSBSliceSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBPhaseSpoilerTable[n] = 0;
                gradAmpSBReadSpoilerTable[n] = grad_amp_sat_spoiler;
                double off_center_neg = parent.getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_2D) - parent.getDouble(CommonUP.FIELD_OF_VIEW_PHASE) / 2.0 + satband_distance_from_fov + satband_thickness / 2.0;
                offsetFreqSBTable[n] = new RFPulse().calculateOffsetFreq(grad_amp_satband_mTpm, off_center_neg);
                n += 1;
            }
            if (position_sli_ph_rea[4] == 1) {
                gradAmpSBSliceTable[n] = 0;
                gradAmpSBPhaseTable[n] = 0;
                gradAmpSBReadTable[n] = grad_amp_satband;
                gradAmpSBSliceSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBPhaseSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBReadSpoilerTable[n] = 0;
                double off_center_pos = parent.getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_1D) + parent.getDouble(CommonUP.FIELD_OF_VIEW) / 2.0 + satband_distance_from_fov + satband_thickness / 2.0;
                offsetFreqSBTable[n] = new RFPulse().calculateOffsetFreq(grad_amp_satband_mTpm, off_center_pos);
                n += 1;
            }
            if (position_sli_ph_rea[5] == 1) {
                gradAmpSBSliceTable[n] = 0;
                gradAmpSBPhaseTable[n] = 0;
                gradAmpSBReadTable[n] = grad_amp_satband;
                gradAmpSBSliceSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBPhaseSpoilerTable[n] = grad_amp_sat_spoiler;
                gradAmpSBReadSpoilerTable[n] = 0;
                double off_center_neg = parent.getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_1D) - parent.getDouble(CommonUP.FIELD_OF_VIEW) / 2.0 + satband_distance_from_fov + satband_thickness / 2.0;
                offsetFreqSBTable[n] = new RFPulse().calculateOffsetFreq(grad_amp_satband_mTpm, off_center_neg);
                n += 1;
            }
        } else if (isTofBandEnabled) {
                double satband_distance_from_slice = parent.getDouble(TofSat.UP.TOF2D_SB_DISTANCE_FROM_SLICE);
                double off_center_slice_pos = satband_distance_from_slice + satband_tof_thickness / 2.0; // sat band cranial from voxel
                double off_center_slice_neg = -off_center_slice_pos;  // caudal
                double off_center_slice = 0;
                if ("BELOW THE SLICE".equalsIgnoreCase(parent.getText(TofSat.UP.TOF2D_SB_POSITION))) {
                    off_center_slice = off_center_slice_neg;
                } else if ("ABOVE THE SLICE".equalsIgnoreCase(parent.getText(TofSat.UP.TOF2D_SB_POSITION))) {
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

    protected void prepPulseComp(boolean isEnabled) {
        pulseTXSatBand.addFrequencyOffset(offsetFreqSBTable);
        pulseTXSatBand.setFrequencyOffset(isEnabled ? Order.LoopB : Order.FourLoop);

        RFPulse pulseTXSatBandPrep = RFPulse.createRFPulse(parent.getSequence(), SP.Time_grad_ramp_sb, SP.Freq_offset_tx_sb_prep);
        pulseTXSatBandPrep.setCompensationFrequencyOffset(pulseTXSatBand, 0.5);
        RFPulse pulseTXSatBandComp = RFPulse.createRFPulse(parent.getSequence(), SP.Time_grad_ramp_sb, SP.Freq_offset_tx_sb_comp);
        pulseTXSatBandComp.setCompensationFrequencyOffset(pulseTXSatBand, 0.5);
    }

    protected void setSeqParamTime() {
        if (isSBorTSEnabled) {
            parent.set(SP.Time_grad_ramp_sb, parent.getDouble(CommonUP.GRADIENT_RISE_TIME));
            parent.set(SP.Time_grad_sb, 0.0005);
            if (parent.hasParam(CommonUP.TX_LENGTH_90)) {
                parent.set(SP.Time_tx_sb, parent.getDouble(CommonUP.TX_LENGTH_90));
            } else if (parent.hasParam(CommonUP.TX_LENGTH)) {
                parent.set(SP.Time_tx_sb, parent.getDouble(CommonUP.TX_LENGTH));
            } else {
                Log.error(getClass(), "User Param TX_LENGTH_90 or TX_LENGTH does not exist");
            }
        } else {
            parent.set(SP.Time_grad_ramp_sb, parent.minInstructionDelay);
            parent.set(SP.Time_grad_sb, parent.minInstructionDelay);
            parent.set(SP.Time_tx_sb, parent.minInstructionDelay);
        }
    }

    protected void clcPulse() {
        if (!pulseTXSatBand.checkPower(getFlipAngle(), parent.observeFrequency, parent.nucleus)) {
//                double tx_length_sb = pulseTXSatBand.getPulseDuration();
//                notifyOutOfRangeParam(TX_LENGTH, pulseTXSatBand.getPulseDuration(), ((NumberParam) getParam(TX_LENGTH)).getMaxValue(), "Pulse length too short to reach RF power with this pulse shape");
            parent.set(SP.Time_tx_sb, pulseTXSatBand.getPulseDuration());
        }
    }

    protected double getFlipAngle() {
        double flip_angle = pulseTXSatBand.isSlr() ? 90 : parent.getDouble(CommonUP.FLIP_ANGLE);
        parent.getParam(CommonUP.FLIP_ANGLE).setValue(flip_angle);

        double flip_angle_satband = 0;
        double time_tau_sat = 0.0;

        if (isSBorTSEnabled) {
            if (parent.hasParam(UP.SATBAND_TAU) && parent.getDouble(UP.SATBAND_TAU) > 0.0) {
                time_tau_sat = parent.getDouble(UP.SATBAND_TAU);
            } else {
                //time_tau_sat = TimeEvents.getTimeBetweenEvents(parent.getSequence(), Events.FatSatPulse.ID, Events.TX90.ID);
                Log.warning(getClass(), "We suppose the time between FatSat and Excitation Pules are negligible");

                if (parent.hasParam(FatSat.UP.FATSAT_TX_LENGTH) && parent.getDouble(FatSat.UP.FATSAT_TX_LENGTH) > 0.0)
                    time_tau_sat = parent.getDouble(FatSat.UP.FATSAT_TX_LENGTH);
                if (parent.hasParam(CommonUP.TX_LENGTH_90))
                    time_tau_sat += parent.getDouble(CommonUP.TX_LENGTH_90);
                if (parent.hasParam(CommonUP.TX_LENGTH))
                    time_tau_sat += parent.getDouble(CommonUP.TX_LENGTH);
            }

            double time_t1_satband = parent.getDouble(UP.SATBAND_T1);
            double t1_relax_time_sat = time_t1_satband / 1000.0;   // T1_tissue = 500ms
            //
            double flip_90_sat = flip_angle == 90 ? Math.acos((1 - Math.exp(time_tau_sat / t1_relax_time_sat)) / (1 - Math.exp((time_tau_sat - parent.getDouble(CommonUP.REPETITION_TIME)) / t1_relax_time_sat))) : Math.acos(1 - Math.exp(time_tau_sat / t1_relax_time_sat));
            double flip_90_sat_degree = Math.toDegrees((isTofBandEnabled ? 1.5 : 1) * flip_90_sat);
            flip_angle_satband = pulseTXSatBand.isSlr() ? 90 : flip_90_sat_degree;  //ha slr,akkor legyen 90,különben szar a szeletprofil!
        }
        return flip_angle_satband;
    }

    protected int[] satBandPrep() {
        String satbandOrientation = parent.getText(UP.SATBAND_ORIENTATION);
        String orientation = parent.getText(CommonUP.ORIENTATION);

        boolean cranial = false;
        boolean caudal = false;
        boolean anterior = false;
        boolean posterior = false;
        boolean right = false;
        boolean left = false;
        if ("CRANIAL".equalsIgnoreCase(satbandOrientation)) {
            cranial = true;
        } else if ("CAUDAL".equalsIgnoreCase(satbandOrientation)) {
            caudal = true;
        } else if ("CRANIAL AND CAUDAL".equalsIgnoreCase(satbandOrientation)) {
            cranial = true;
            caudal = true;
        } else if ("ANTERIOR".equalsIgnoreCase(satbandOrientation)) {
            anterior = true;
        } else if ("POSTERIOR".equalsIgnoreCase(satbandOrientation)) {
            posterior = true;
        } else if ("ANTERIOR AND POSTERIOR".equalsIgnoreCase(satbandOrientation)) {
            anterior = true;
            posterior = true;
        } else if ("RIGHT".equalsIgnoreCase(satbandOrientation)) {
            right = true;
        } else if ("LEFT".equalsIgnoreCase(satbandOrientation)) {
            left = true;
        } else if ("RIGHT AND LEFT".equalsIgnoreCase(satbandOrientation)) {
            right = true;
            left = true;
        } else if ("RIGHT AND LEFT".equalsIgnoreCase(satbandOrientation)) {
            right = true;
            left = true;
        } else if ("ALL".equalsIgnoreCase(satbandOrientation)) {
            cranial = true;
            caudal = true;
            anterior = true;
            posterior = true;
            right = true;
            left = true;
        }

        position_sli_ph_rea[0] = 0;
        position_sli_ph_rea[1] = 0;
        position_sli_ph_rea[2] = 0;
        position_sli_ph_rea[3] = 0;
        position_sli_ph_rea[4] = 0;
        position_sli_ph_rea[5] = 0;
        if ("AXIAL".equalsIgnoreCase(orientation)) {
            position_sli_ph_rea[0] = cranial ? 1 : 0;
            position_sli_ph_rea[1] = caudal ? 1 : 0;
            position_sli_ph_rea[2] = anterior ? 1 : 0;
            position_sli_ph_rea[3] = posterior ? 1 : 0;
            position_sli_ph_rea[4] = right ? 1 : 0;
            position_sli_ph_rea[5] = left ? 1 : 0;
        } else if ("SAGITTAL".equalsIgnoreCase(orientation)) {
            position_sli_ph_rea[0] = left ? 1 : 0;
            position_sli_ph_rea[1] = right ? 1 : 0;
            position_sli_ph_rea[2] = anterior ? 1 : 0;
            position_sli_ph_rea[3] = posterior ? 1 : 0;
            position_sli_ph_rea[4] = cranial ? 1 : 0;
            position_sli_ph_rea[5] = caudal ? 1 : 0;
        } else if ("CORONAL".equalsIgnoreCase(orientation)) {
            position_sli_ph_rea[0] = anterior ? 1 : 0;
            position_sli_ph_rea[1] = posterior ? 1 : 0;
            position_sli_ph_rea[2] = right ? 1 : 0;
            position_sli_ph_rea[3] = left ? 1 : 0;
            position_sli_ph_rea[4] = cranial ? 1 : 0;
            position_sli_ph_rea[5] = caudal ? 1 : 0;
        } else if ("OBLIQUE".equalsIgnoreCase(orientation)) {
            List<Double> image_orientation = parent.getListDouble(CommonUP.IMAGE_ORIENTATION_SUBJECT);
            double[][] dir_ind = new double[3][3];
            for (int i = 0; i < 3; i++) {
                dir_ind[0][i] = image_orientation.get(i);
                dir_ind[1][i] = image_orientation.get(i + 3);
            }
            dir_ind[2][0] = dir_ind[0][1] * dir_ind[1][2] - dir_ind[0][2] * dir_ind[1][1];
            dir_ind[2][1] = dir_ind[0][2] * dir_ind[1][0] - dir_ind[0][0] * dir_ind[1][2];
            dir_ind[2][2] = dir_ind[0][0] * dir_ind[1][1] - dir_ind[0][1] * dir_ind[1][0];
            int i, j;
            int max_index = 0;
            double norm_vector_re = Math.sqrt(Math.pow(dir_ind[0][0], 2) + Math.pow(dir_ind[0][1], 2) + Math.pow(dir_ind[0][2], 2));
            double norm_vector_ph = Math.sqrt(Math.pow(dir_ind[1][0], 2) + Math.pow(dir_ind[1][1], 2) + Math.pow(dir_ind[1][2], 2));
            double norm_vector_sl = Math.sqrt(Math.pow(dir_ind[2][0], 2) + Math.pow(dir_ind[2][1], 2) + Math.pow(dir_ind[2][2], 2));
            //normalizing vectors
            dir_ind[0][0] = dir_ind[0][0] / norm_vector_re;
            dir_ind[0][1] = dir_ind[0][1] / norm_vector_re;
            dir_ind[0][2] = dir_ind[0][2] / norm_vector_re;
            dir_ind[1][0] = dir_ind[1][0] / norm_vector_ph;
            dir_ind[1][1] = dir_ind[1][1] / norm_vector_ph;
            dir_ind[1][2] = dir_ind[1][2] / norm_vector_ph;
            dir_ind[2][0] = dir_ind[2][0] / norm_vector_sl;
            dir_ind[2][1] = dir_ind[2][1] / norm_vector_sl;
            dir_ind[2][2] = dir_ind[2][2] / norm_vector_sl;
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    System.out.println("dir_ind[" + i + "][" + j + "]" + dir_ind[i][j]);
                }
            }

            // System.out.println(" direction index and dir ind:  "+direction_index[2]+" "+dir_ind[0][2]);
            int[] max_vector = new int[3];

            // read, phase and slice vector which component has the largest value
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    if (Math.abs(dir_ind[i][j]) >= Math.abs(dir_ind[i][max_index])) {
                        max_index = j;
                    }
                }
                max_vector[i] = max_index; // storing each vector's maximum value index
                //System.out.println("max_vector["+i+"]"+max_vector[i]);
            }

            boolean[][] anatomy_to_local_mx = new boolean[6][6];

            for (i = 0; i < 6; i++) {
                for (j = 0; j < 6; j++) {
                    anatomy_to_local_mx[i][j] = false;
                }
            }

            for (i = 0; i < 3; i++) {

                if (dir_ind[i][max_vector[i]] < 0) {
                    anatomy_to_local_mx[i][max_vector[i] + 3] = true;
                    anatomy_to_local_mx[i + 3][max_vector[i]] = true;
                } else {
                    anatomy_to_local_mx[i][max_vector[i]] = true;
                    anatomy_to_local_mx[i + 3][max_vector[i] + 3] = true;
                }
            }
            boolean[] local_vector = new boolean[6];

            local_vector[0] = false;
            local_vector[1] = false;
            local_vector[2] = false;
            local_vector[3] = false;
            local_vector[4] = false;
            local_vector[5] = false;

            boolean[] anatomy_vector = new boolean[6];

            anatomy_vector[0] = right;
            anatomy_vector[1] = posterior;
            anatomy_vector[2] = caudal;
            anatomy_vector[3] = left;
            anatomy_vector[4] = anterior;
            anatomy_vector[5] = cranial;

            boolean sum;
            for (i = 0; i < 6; i++) {
                sum = false;
                for (j = 0; j < 6; j++) {
                    sum = sum || (anatomy_to_local_mx[i][j] & anatomy_vector[j]);
                    //	System.out.println("sum= "+sum+" + "+anatomy_to_local_mx[i][j]+"*"+anatomy_vector[j]);
                }
                local_vector[i] = sum;
                // System.out.println("local vector "+local_vector[i]);

            }
            position_sli_ph_rea[4] = local_vector[0] ? 1 : 0;
            position_sli_ph_rea[2] = local_vector[1] ? 1 : 0;
            position_sli_ph_rea[0] = local_vector[2] ? 1 : 0;
            position_sli_ph_rea[5] = local_vector[3] ? 1 : 0;
            position_sli_ph_rea[3] = local_vector[4] ? 1 : 0;
            position_sli_ph_rea[1] = local_vector[5] ? 1 : 0;

            // System.out.println("read+ "+position_sli_ph_rea[4]+" phase+ "+position_sli_ph_rea[2]+" slice+ "+position_sli_ph_rea[0]);
            // System.out.println("read- "+position_sli_ph_rea[5]+" phase- "+position_sli_ph_rea[3]+" slice- "+position_sli_ph_rea[1]);
        }
        boolean is_switch = parent.getBoolean(CommonUP.SWITCH_READ_PHASE);
        boolean phase_pos_temp = position_sli_ph_rea[2] == 1;
        boolean phase_neg_temp = position_sli_ph_rea[3] == 1;
        boolean read_pos_temp = position_sli_ph_rea[4] == 1;
        boolean read_neg_temp = position_sli_ph_rea[5] == 1;
        if (is_switch) {
            position_sli_ph_rea[2] = read_pos_temp ? 1 : 0;
            position_sli_ph_rea[3] = read_neg_temp ? 1 : 0;
            position_sli_ph_rea[4] = phase_pos_temp ? 1 : 0;
            position_sli_ph_rea[5] = phase_neg_temp ? 1 : 0;
        }
        return position_sli_ph_rea;
    }

}

