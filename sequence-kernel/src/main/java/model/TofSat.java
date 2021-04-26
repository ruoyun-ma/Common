package model;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;

import common.*;
import kernel.*;

import java.util.ArrayList;
import java.util.List;

import static common.CommonUP.*;
import static common.CommonSP.*;

/**
 * The Current Design does not allow one to use TofSat and SatBand
 * at same time, But one may improve it afterwards (if necessary)
 */


public class TofSat extends SatBand {

    protected enum UP implements GeneratorParamEnum {
        TOF2D_ENABLED,
        TOF2D_SB_THICKNESS,
        TOF2D_SB_TX_SHAPE,
        TOF2D_SB_GRAMP_SP,
        TOF2D_SB_DISTANCE_FROM_SLICE,
        TOF2D_SB_POSITION,
        TOF2D_SB_OFFSET,
        TOF3D_MT_FLIP_ANGLE,
        TOF3D_TX_RAMP_SLOPE,
        TOF3D_MOTSA_OVERLAP,
        TOF3D_EXT_SHIRNK_FACTOR,
        TOF3D_MT_TX_LENGTH,
        TOF2D_FLOW_VELOCITY,
        TOF2D_FLOW_TAU,
        ;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Time_flow,
        ;
    }

    public TofSat(SeqPrep parent) {
        super(parent);
        this.parent = parent;
    }

    @Override
    public void init() {
        super.init();
        //isTofBandEnabled = parent.getBoolean(MULTI_PLANAR_EXCITATION) && parent.getBoolean(UP.TOF2D_ENABLED); // TOF not allowed in 3D
        isTofBandEnabled = parent.getBoolean(UP.TOF2D_ENABLED); //XG: now we do support 3D

        if (isTofBandEnabled) { //TofBandEnabled has a higher priority than SatBandEnabled
            isSatBandEnabled = false;
            if (parent.hasParam(SatBand.UP.SATBAND_ENABLED)) {
                parent.getParam(SatBand.UP.SATBAND_ENABLED).setValue(isSatBandEnabled);
            }
        }
        List<String> tmp_tx_shape = new ArrayList<>(parent.tx_shape);
        tmp_tx_shape.add("RAMP");
        parent.tx_shape = tmp_tx_shape;
    }

    @Override
    protected void prepGrad() {
        super.prepGrad();
        getFlow();
    }

    @Override
    public String getName() {
        return "TofSat";
    }

    @Override
    public boolean isEnabled() {
        return isTofBandEnabled;
    }

    private void getFlow() {
        double sat_flow_time_corr = parent.minInstructionDelay;
        if (isTofBandEnabled) {
            double satband_distance_from_slice = parent.getDouble(UP.TOF2D_SB_DISTANCE_FROM_SLICE);
            double sat_flow_dist = satband_distance_from_slice + (parent.isMultiplanar ? parent.sliceThickness : (parent.sliceThickness * parent.userMatrixDimension3D));
            double sat_flow_velocity = parent.getDouble(UP.TOF2D_FLOW_VELOCITY);
            double sat_flow_time = sat_flow_dist / sat_flow_velocity;
            double time = 0.0;
            if (parent.hasParam(UP.TOF2D_FLOW_TAU) && parent.roundToDecimal(parent.getDouble(UP.TOF2D_FLOW_TAU),6) >= parent.minInstructionDelay) {
                time = parent.getDouble(UP.TOF2D_FLOW_TAU);
                sat_flow_time_corr = parent.getDouble(UP.TOF2D_FLOW_TAU);
            } else {
                //                TimeEvents.getTimeBetweenEvents(parent.getSequence(), Events.LoopSatBandEnd.ID - 1, Events.P90.ID - 3)
//                        + TimeEvents.getTimeBetweenEvents(getSequence(), Events.P90.ID - 2, Events.Acq.ID);//(real index - 1) in the argument instead of real index
                time += parent.getDouble(GRADIENT_RISE_TIME)
                        + (parent.hasParam(TX_LENGTH_90) ? parent.getDouble(TX_LENGTH_90) : parent.getDouble(TX_LENGTH)) / 2
                        + parent.getDouble(ECHO_TIME)
                        + parent.getSequenceTable(Time_rx).get(0).doubleValue()/2;

                if (parent.models.containsKey("SatBand")) {
                    time += parent.getSequenceTable(SatBand.SP.Time_delay_sb).get(0).doubleValue();
                }
                if (parent.models.containsKey("FatSat") && parent.models.containsKey("FatSatWep")) {
                    time += parent.models.get("FatSatWep").getDuration();
                } else if (parent.models.containsKey("FatSat")) {
                    time += parent.models.get("FatSat").getDuration();
                }
                if (parent.models.containsKey("InvRec")) {
                    time += parent.getSequenceTable(InvRec.SP.Time_TI_delay).getMaxValue();
                }
//                System.out.println("time" + time);
//                System.out.println("sat_flow_time" + sat_flow_time);
                sat_flow_time_corr = sat_flow_time - time;
                sat_flow_time_corr = (sat_flow_time_corr < 0) ? 0.0001 : sat_flow_time_corr;
//                System.out.println("sat_flow_time_corr" + sat_flow_time_corr);
            }
        }
        parent.set(SP.Time_flow, sat_flow_time_corr);
    }

}

