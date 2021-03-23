package rs2d.sequence.common;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;

/**
 *  V1.0 - 2021.3 XG
 *
 */

public class FatSat implements ModelInterface {
    public SeqPrep parent;
    protected static boolean isFatSatEnabled;
    protected static boolean isFatSatWepEnabled;

    public RFPulse pulseTXFatSat;
    public Gradient gradFatsatRead;
    public Gradient gradFatsatPhase;
    public Gradient gradFatsatSlice;

    protected boolean isAttAuto;
    private boolean isFSorFSWEnabled;

    //TODO: I really thought we should unify the name one day
    protected enum UP implements GeneratorParamEnum {
        FAT_SATURATION_ENABLED,
        FATSAT_GRAD_APP_TIME,
        FATSAT_OFFSET_FREQ,
        FATSAT_TX_SHAPE,
        FATSAT_TX_LENGTH,
        FATSAT_BANDWIDTH,
        FATSAT_TX_AMP_90,
        FATSAT_TX_AMP,
        FATSAT_FLIP_ANGLE,
        ;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Enable_fs,
        Grad_amp_fatsat_read,
        Grad_amp_fatsat_phase,
        Grad_amp_fatsat_slice,
        Time_grad_fatsat,
        Time_grad_ramp_fatsat,
        Time_before_fatsat_pulse,
        Freq_offset_tx_fatsat_prep,
        Freq_offset_tx_fatsat_comp,
        Tx_amp_fatsat,
        Tx_phase_fatsat,
        Time_tx_fatsat,
        Tx_shape_fatsat,
        Tx_shape_phase_fatsat,
        Freq_offset_tx_fatsat,
        ;
    }

    public FatSat(SeqPrep parent) {
        this.parent = parent;
    }

    @Override
    public void init() {
        parent.setSuggestedValFromListString(parent.tx_shape, true, UP.FATSAT_TX_SHAPE);

        isAttAuto = parent.getBoolean(CommonUP.TX_AMP_ATT_AUTO);
        isFatSatEnabled = parent.getBoolean(UP.FAT_SATURATION_ENABLED);
    }

    @Override
    public void initFinal() throws Exception {
        isFSorFSWEnabled = isFatSatEnabled || isFatSatWepEnabled;
        parent.set(SP.Enable_fs, isFSorFSWEnabled);

        parent.getParam(UP.FATSAT_BANDWIDTH).setDefaultValue(parent.protonFrequency * 3.5);
        parent.getParam(UP.FATSAT_OFFSET_FREQ).setDefaultValue(parent.protonFrequency * 3.5);

        if (parent.hasParam(UP.FAT_SATURATION_ENABLED)) {
            parent.getParam(UP.FAT_SATURATION_ENABLED).setValue(isFatSatEnabled);
        }
        if (parent.hasParam(FatSatWep.UP.FAT_SATURATION_WEP_ENABLED)) {
            parent.getParam(FatSatWep.UP.FAT_SATURATION_WEP_ENABLED).setValue(isFatSatWepEnabled);
        }

        setSeqParamTime();
        initPulseandGrad();
    }

    @Override
    public void prep() throws Exception {
        prepPulse();
        prepGrad();
        prepPulseComp();
    }

    @Override
    public void prepFinal() {
        if (isAttAuto) {
            if (parent.hasParam(UP.FATSAT_TX_AMP_90)) {
                parent.getParam(UP.FATSAT_TX_AMP_90).setValue(pulseTXFatSat.getAmp90());
            } else if (parent.hasParam(UP.FATSAT_TX_AMP)) {
                parent.getParam(UP.FATSAT_TX_AMP).setValue(pulseTXFatSat.getAmp90());
            }
        } else {
            if (parent.hasParam(UP.FATSAT_TX_AMP_90)) {
                pulseTXFatSat.setAmp(parent.getDouble(UP.FATSAT_TX_AMP_90));
            } else if (parent.hasParam(UP.FATSAT_TX_AMP)) {
                pulseTXFatSat.setAmp(parent.getDouble(UP.FATSAT_TX_AMP));
            }
        }

        if (parent.hasParam(CommonUP.SEQ_DESCRIPTION)) {
            String seqDescription = parent.getText(CommonUP.SEQ_DESCRIPTION);
            if (isFatSatEnabled) {
                seqDescription += "_FATSAT";
            }
            if (isFatSatWepEnabled) {
                seqDescription += "_FATSATWEP";
            }
            parent.getParam(CommonUP.SEQ_DESCRIPTION).setValue(seqDescription);
        }
    }

    @Override
    public RFPulse getRfPulses() {
        return pulseTXFatSat;
    }

    @Override
    public double getDuration() {
        double GradDuration = 2 * parent.getSequenceTable(SP.Time_grad_ramp_fatsat).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_grad_fatsat).get(0).doubleValue();

        double RFDuration = parent.getSequenceTable(SP.Time_before_fatsat_pulse).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_tx_fatsat).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_grad_fatsat).get(0).doubleValue();

        double OverlapDuration = parent.getSequenceTable(SP.Time_grad_ramp_fatsat).get(0).doubleValue();

        double eachDuration = GradDuration + RFDuration - OverlapDuration;

        return eachDuration;
    }

    @Override
    public String getName() {
        return "FatSat";
    }

    @Override
    public boolean isEnabled() {
        return isFatSatEnabled;
    }

    protected void initPulseandGrad() throws Exception {
        pulseTXFatSat = RFPulse.createRFPulse(parent.getSequence(), CommonSP.Tx_att, SP.Tx_amp_fatsat, SP.Tx_phase_fatsat,
                SP.Time_tx_fatsat, SP.Tx_shape_fatsat, SP.Tx_shape_phase_fatsat, SP.Freq_offset_tx_fatsat);
        pulseTXFatSat.setShape(parent.getText(UP.FATSAT_TX_SHAPE), parent.nb_shape_points, "Hamming");

        gradFatsatRead = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_fatsat_read, SP.Time_grad_fatsat,
                CommonSP.Grad_shape_rise_up, CommonSP.Grad_shape_rise_down, SP.Time_grad_ramp_fatsat, parent.nucleus);
        gradFatsatPhase = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_fatsat_phase, SP.Time_grad_fatsat,
                CommonSP.Grad_shape_rise_up, CommonSP.Grad_shape_rise_down, SP.Time_grad_ramp_fatsat, parent.nucleus);
        gradFatsatSlice = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_fatsat_slice, SP.Time_grad_fatsat,
                CommonSP.Grad_shape_rise_up, CommonSP.Grad_shape_rise_down, SP.Time_grad_ramp_fatsat, parent.nucleus);
    }

    protected void prepPulse() throws Exception {
        // Fat SAT RF pulse
        parent.getParam(UP.FATSAT_FLIP_ANGLE).setValue(getFlipAngle());

        if (isAttAuto)
            clcPulse();
    }

    protected void prepGrad() {
        if (isFSorFSWEnabled) {
            double pixelDimension = parent.getDouble(CommonUP.RESOLUTION_FREQUENCY);
            double pixel_dimension_ph = parent.getDouble(CommonUP.RESOLUTION_PHASE);
            double pixel_dimension_sl = parent.getDouble(CommonUP.RESOLUTION_SLICE);

            int pixmax = (pixelDimension > pixel_dimension_ph) ? 1 : 2;
            pixmax = (Math.max(pixelDimension, pixel_dimension_ph) > pixel_dimension_sl) ? pixmax : 3;
            boolean test_grad;
            double min_fatsat_application_time;
            switch (pixmax) {
                case 1:
                    test_grad = gradFatsatRead.addSpoiler(pixelDimension, 3);
                    min_fatsat_application_time = gradFatsatRead.getMinTopTime();
                    break;
                case 2:
                    test_grad = gradFatsatPhase.addSpoiler(pixel_dimension_ph, 3);
                    min_fatsat_application_time = gradFatsatPhase.getMinTopTime();
                    break;
                default:
                    test_grad = gradFatsatSlice.addSpoiler(pixel_dimension_sl, 3);
                    min_fatsat_application_time = gradFatsatSlice.getMinTopTime();
            }

            if (!test_grad) {
                parent.notifyOutOfRangeParam(UP.FATSAT_GRAD_APP_TIME, min_fatsat_application_time, ((NumberParam) parent.getParam(UP.FATSAT_GRAD_APP_TIME)).getMaxValue(), "FATSAT_GRAD_APP_TIME too short to get correct Spoiling");
                parent.set(SP.Time_grad_fatsat, min_fatsat_application_time);
                gradFatsatRead.rePrepare();
                gradFatsatPhase.rePrepare();
                gradFatsatSlice.rePrepare();
            }
        }

        gradFatsatRead.applyAmplitude();
        gradFatsatPhase.applyAmplitude();
        gradFatsatSlice.applyAmplitude();
    }

    protected void prepPulseComp() {
        pulseTXFatSat.setFrequencyOffset(isFatSatEnabled ? parent.getDouble(UP.FATSAT_OFFSET_FREQ) : 0.0);

        RFPulse pulseTXFatSatPrep = RFPulse.createRFPulse(parent.getSequence(), SP.Time_before_fatsat_pulse, SP.Freq_offset_tx_fatsat_prep);
        pulseTXFatSatPrep.setCompensationFrequencyOffset(pulseTXFatSat, 0.5);
        RFPulse pulseTXFatSatComp = RFPulse.createRFPulse(parent.getSequence(), SP.Time_grad_ramp_fatsat, SP.Freq_offset_tx_fatsat_comp);
        pulseTXFatSatComp.setCompensationFrequencyOffset(pulseTXFatSat, 0.5);
    }

    protected void setSeqParamTime() {
        double tx_bandwidth_90_fs = parent.getDouble(UP.FATSAT_BANDWIDTH);
        double tx_bandwidth_factor_90_fs = parent.getTx_bandwidth_factor(UP.FATSAT_TX_SHAPE, CommonUP.TX_BANDWIDTH_FACTOR, CommonUP.TX_BANDWIDTH_FACTOR_3D);
        double tx_length_90_fs = isFatSatEnabled ? tx_bandwidth_factor_90_fs / tx_bandwidth_90_fs : parent.minInstructionDelay;
        tx_length_90_fs = isFatSatWepEnabled ? parent.getDouble(FatSatWep.UP.FATSAT_WEP_TX_LENGTH) : tx_length_90_fs;

        parent.getParam(UP.FATSAT_TX_LENGTH).setValue(tx_length_90_fs);
        parent.set(SP.Time_tx_fatsat, tx_length_90_fs);
        parent.set(SP.Time_before_fatsat_pulse, parent.blankingDelay);

        // FatSAT timing seting needed for Flip Angle calculation
        if (isFSorFSWEnabled) {
            parent.set(SP.Time_grad_fatsat, parent.getDouble(UP.FATSAT_GRAD_APP_TIME));
            parent.set(SP.Time_grad_ramp_fatsat, parent.getDouble(CommonUP.GRADIENT_RISE_TIME));
        } else {
            parent.set(SP.Time_grad_fatsat, parent.minInstructionDelay);
            parent.set(SP.Time_grad_ramp_fatsat, parent.minInstructionDelay);
        }

        if (isFatSatWepEnabled) {
            parent.set(SP.Time_tx_fatsat, parent.getDouble(FatSatWep.UP.FATSAT_WEP_TX_LENGTH));
        }
    }

    protected void clcPulse() {
        if (!pulseTXFatSat.checkPower(getFlipAngle(), getFrequency(), parent.nucleus)) {
            double tx_length_90_fs = pulseTXFatSat.getPulseDuration();
            parent.set(SP.Time_tx_fatsat, tx_length_90_fs);
            parent.getParam(UP.FATSAT_TX_LENGTH).setValue(pulseTXFatSat.getPulseDuration());
            if (isFatSatWepEnabled) {
                //parent.set(SP.Time_tx_fatsat, tx_length_90_fs);
                parent.set(FatSatWep.SP.Time_tx_fatsat_wep, tx_length_90_fs);
                parent.getParam(FatSatWep.UP.FATSAT_WEP_TX_LENGTH).setValue(tx_length_90_fs);
            }
            System.out.println(pulseTXFatSat.getPower()+ " pulseTXFatSat.getPower()   ERROR");
        }
    }

    protected double getFlipAngle() {
        if (isFatSatEnabled)
            return 90.0;
        else if (isFatSatWepEnabled)
            return 45.0;
        else
            return 0.0;
    }

    protected double getFrequency() {
        if (isFatSatEnabled)
            return parent.observeFrequency + parent.getDouble(UP.FATSAT_OFFSET_FREQ);
        else if (isFatSatWepEnabled)
            return parent.observeFrequency;
        else
            return parent.observeFrequency;
    }

}
