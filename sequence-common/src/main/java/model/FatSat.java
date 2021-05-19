package model;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;

import common.*;
import kernel.*;

import static common.CommonUP.*;
import static common.CommonSP.*;

/**
 * V1.0 - 2021.3 XG
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
        FATSAT_DELAY,
        ;

        @Override
        public Param build() {
            throw new UnsupportedOperationException("Unable to create a Param");
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Enable_fs,
        Tx_att_offset_fs,
        Tx_amp_fs,
        Tx_phase_fs,
        Tx_shape_fs,
        Tx_shape_phase_fs,
        Time_tx_fs,
        Time_grad_fs,
        Time_grad_ramp_fs,
        Time_delay_fs,
        Time_before_fs_pulse,
        Freq_offset_tx_fs,
        Freq_offset_tx_fs_prep,
        Freq_offset_tx_fs_comp,

        Grad_amp_fs_read,
        Grad_amp_fs_phase,
        Grad_amp_fs_slice,
        ;
    }

    public FatSat() {
    }

    public FatSat(SeqPrep parent) {
        this.parent = parent;
    }

    public void init() {
        parent.setSuggestedValFromListString(parent.tx_shape, true, UP.FATSAT_TX_SHAPE);

        isAttAuto = parent.getBoolean(TX_AMP_ATT_AUTO);
        isFatSatEnabled = parent.getBoolean(UP.FAT_SATURATION_ENABLED);
    }

    @Override
    public void init(SeqPrep parent) {
        this.parent = parent;
        init();
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

        if (parent.hasParam(SEQ_DESCRIPTION)) {
            String seqDescription = parent.getText(SEQ_DESCRIPTION);
            if (isFatSatEnabled) {
                seqDescription += "_FATSAT";
            }
            if (isFatSatWepEnabled) {
                seqDescription += "_FATSATWEP";
            }
            parent.getParam(SEQ_DESCRIPTION).setValue(seqDescription);
        }
    }

    @Override
    public RFPulse getRfPulses() {
        return pulseTXFatSat;
    }

    @Override
    public double getDuration() {
        double GradDuration = 2 * parent.getSequenceTable(SP.Time_grad_ramp_fs).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_grad_fs).get(0).doubleValue();
        double RFDuration = parent.getSequenceTable(SP.Time_before_fs_pulse).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_tx_fs).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_grad_fs).get(0).doubleValue();
        double Delay = 0.0;
        if (parent.getSequence().getPublicTable(SP.Time_delay_fs.name()) != null)
            Delay = parent.getSequenceTable(SP.Time_delay_fs).get(0).doubleValue();

        double OverlapDuration = parent.getSequenceTable(SP.Time_grad_ramp_fs).get(0).doubleValue();
        double eachDuration = GradDuration + RFDuration - OverlapDuration + Delay;

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
        pulseTXFatSat = RFPulse.createRFPulse(parent.getSequence(), Tx_att, SP.Tx_amp_fs, SP.Tx_phase_fs,
                SP.Time_tx_fs, SP.Tx_shape_fs, SP.Tx_shape_phase_fs, SP.Freq_offset_tx_fs);
        if (parent.getSequence().getPublicTable(SP.Tx_att_offset_fs.name()) != null) {
            pulseTXFatSat.createAttOffset(parent.getSequence(), SP.Tx_att_offset_fs);
        }

        pulseTXFatSat.setShape(parent.getText(UP.FATSAT_TX_SHAPE), parent.nb_shape_points, "Hamming");

        gradFatsatRead = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_fs_read, SP.Time_grad_fs,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_fs, parent.nucleus);
        gradFatsatPhase = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_fs_phase, SP.Time_grad_fs,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_fs, parent.nucleus);
        gradFatsatSlice = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_fs_slice, SP.Time_grad_fs,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_fs, parent.nucleus);
    }

    protected void prepPulse() throws Exception {
        // Fat SAT RF pulse
        parent.getParam(UP.FATSAT_FLIP_ANGLE).setValue(getFlipAngle());

        if (isAttAuto)
            clcPulse();
    }

    protected void prepGrad() {
        if (isFSorFSWEnabled) {
            double pixelDimension = parent.getDouble(RESOLUTION_FREQUENCY);
            double pixel_dimension_ph = parent.getDouble(RESOLUTION_PHASE);
            double pixel_dimension_sl = parent.getDouble(RESOLUTION_SLICE);

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
                parent.set(SP.Time_grad_fs, min_fatsat_application_time);
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

        RFPulse pulseTXFatSatPrep = RFPulse.createRFPulse(parent.getSequence(), SP.Time_before_fs_pulse, SP.Freq_offset_tx_fs_prep);
        pulseTXFatSatPrep.setCompensationFrequencyOffset(pulseTXFatSat, 0.5);
        RFPulse pulseTXFatSatComp = RFPulse.createRFPulse(parent.getSequence(), SP.Time_grad_ramp_fs, SP.Freq_offset_tx_fs_comp);
        pulseTXFatSatComp.setCompensationFrequencyOffset(pulseTXFatSat, 0.5);
    }

    protected void setSeqParamTime() {
        if (parent.getSequence().getPublicTable(SP.Time_delay_fs.name()) != null){
        if (parent.hasParam(UP.FATSAT_DELAY))
            parent.set(SP.Time_delay_fs, UP.FATSAT_DELAY);
        else
            parent.set(SP.Time_delay_fs, parent.minInstructionDelay);
        }

        double tx_bandwidth_90_fs = parent.getDouble(UP.FATSAT_BANDWIDTH);
        double tx_bandwidth_factor_90_fs = parent.getTx_bandwidth_factor(UP.FATSAT_TX_SHAPE, TX_BANDWIDTH_FACTOR, TX_BANDWIDTH_FACTOR_3D);
        double tx_length_90_fs = isFatSatEnabled ? tx_bandwidth_factor_90_fs / tx_bandwidth_90_fs : parent.minInstructionDelay;
        tx_length_90_fs = isFatSatWepEnabled ? parent.getDouble(FatSatWep.UP.FATSAT_WEP_TX_LENGTH) : tx_length_90_fs;

        parent.getParam(UP.FATSAT_TX_LENGTH).setValue(tx_length_90_fs);
        parent.set(SP.Time_tx_fs, tx_length_90_fs);
        parent.set(SP.Time_before_fs_pulse, parent.blankingDelay);

        // FatSAT timing seting needed for Flip Angle calculation
        if (isFSorFSWEnabled) {
            parent.set(SP.Time_grad_fs, parent.getDouble(UP.FATSAT_GRAD_APP_TIME));
            parent.set(SP.Time_grad_ramp_fs, parent.getDouble(GRADIENT_RISE_TIME));
        } else {
            parent.set(SP.Time_grad_fs, parent.minInstructionDelay);
            parent.set(SP.Time_grad_ramp_fs, parent.minInstructionDelay);
        }

        if (isFatSatWepEnabled) {
            parent.set(SP.Time_tx_fs, parent.getDouble(FatSatWep.UP.FATSAT_WEP_TX_LENGTH));
        }
    }

    protected void clcPulse() {
        if (pulseTXFatSat == null){
            System.out.println("00000000000000");

        }
        if (!pulseTXFatSat.checkPower(getFlipAngle(), getFrequency(), parent.nucleus)) {
            double tx_length_90_fs = pulseTXFatSat.getPulseDuration();
            parent.set(SP.Time_tx_fs, tx_length_90_fs);
            parent.getParam(UP.FATSAT_TX_LENGTH).setValue(pulseTXFatSat.getPulseDuration());
            if (isFatSatWepEnabled) {
                //parent.set(SP.Time_tx_fatsat, tx_length_90_fs);
                parent.set(FatSatWep.SP.Time_tx_fs_wep, tx_length_90_fs);
                parent.getParam(FatSatWep.UP.FATSAT_WEP_TX_LENGTH).setValue(tx_length_90_fs);
            }
            System.out.println(pulseTXFatSat.getPower() + " pulseTXFatSat.getPower()   ERROR");
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
