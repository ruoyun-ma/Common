package model;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;

import kernel.*;

/**
 * V1.0 - 2021.3 XG
 */

public class FatSatWep extends FatSat {

    //TODO: I really thought we should unify the name one day
    protected enum UP implements GeneratorParamEnum {
        FAT_SATURATION_WEP_ENABLED,
        FATSAT_WEP_TX_LENGTH,
        ;

        @Override
        public Param build() {
            throw new UnsupportedOperationException("Unable to create a Param");
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Enable_fs_wep,
        Time_tx_fs_wep,
        Time_before_fs_wep_pulse,
        Time_delay_wep,
        Tx_phase_fs_wep,
        Freq_offset_tx_fs_wep,
        ;
    }

    public FatSatWep() {
    }

    @Override
    public void init(SeqPrep parent) {
        super.init(parent);
    }

    @Override
    public void initPre() throws Exception {
        super.initPre();
        isFatSatWepEnabled = parent.getBoolean(UP.FAT_SATURATION_WEP_ENABLED);

        if (isFatSatWepEnabled) { //FatSatWepEnabled has a higher priority than FatSatEnabled
            isFatSatEnabled = false;
            if (parent.hasParam(FatSat.UP.FAT_SATURATION_ENABLED)) {
                parent.getParam(FatSat.UP.FAT_SATURATION_ENABLED).setValue(isFatSatEnabled);
            }
        }
    }

    @Override
    public double getDuration() {
        double FatSat = super.getDuration();

        double RFFatWepDuration = parent.getSequenceTable(SP.Time_before_fs_wep_pulse).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_tx_fs_wep).get(0).doubleValue();

        double FatWepDelayDuration = parent.getSequenceTable(SP.Time_delay_wep).get(0).doubleValue();

        return FatSat + RFFatWepDuration + FatWepDelayDuration;
    }

    @Override
    public void initFinal() throws Exception {
        super.initFinal();
        parent.set(SP.Enable_fs_wep, isFatSatWepEnabled);
    }

    @Override
    public String getName() {
        return "FatSatWep";
    }

    @Override
    public boolean isEnabled() {
        return isFatSatWepEnabled;
    }

    @Override
    protected void setSeqParamTime() {
        super.setSeqParamTime();

        parent.set(SP.Time_before_fs_wep_pulse, parent.blankingDelay);

        if (isFatSatWepEnabled) {
            parent.set(SP.Time_tx_fs_wep, parent.getDouble(UP.FATSAT_WEP_TX_LENGTH));
        } else {
            parent.set(SP.Time_tx_fs_wep, parent.minInstructionDelay);
        }

        double tx_frequency_offset_90_fs = parent.getDouble(FatSat.UP.FATSAT_OFFSET_FREQ);
        double timeFatSatWepDelay = parent.minInstructionDelay;
        double phaseFatSat = 0.0;
        if (isFatSatWepEnabled) {
            double timeFatSatWepBtwPulse = (parent.getSequenceTable(FatSat.SP.Time_tx_fs).get(0).doubleValue()
                    + parent.getSequenceTable(SP.Time_tx_fs_wep).get(0).doubleValue()) / 2
                    + parent.getSequenceTable(FatSat.SP.Time_before_fs_pulse).get(0).doubleValue();

            int nbFatWaterRotation = (int) Math.round(timeFatSatWepBtwPulse / (1 / Math.abs(tx_frequency_offset_90_fs)));
            timeFatSatWepDelay = (1 / Math.abs(tx_frequency_offset_90_fs)) * (nbFatWaterRotation + 1 / 2.0)
                    - timeFatSatWepBtwPulse;

            phaseFatSat = 180;
        }

        parent.set(SP.Time_delay_wep, timeFatSatWepDelay);
        parent.set(SP.Tx_phase_fs_wep, 0);
        parent.set(FatSat.SP.Tx_phase_fs, phaseFatSat);
    }
}
