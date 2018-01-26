package rs2d.sequence.common;

import rs2d.spinlab.instrument.Instrument;

//import rs2d.sequence.common.HardwareShim;
//         HardwareShim hardwareShim = new HardwareShim();
//         setParamValue(HARDWARE_SHIM, hardwareShim.getValue());
//         setParamValue(HARDWARE_SHIM_LABEL, hardwareShim.getLabel());

public class HardwareGradientDelay {
    private double rxTxDelay;
    private double gradientDelay;

    public HardwareGradientDelay() throws Exception {
        rxTxDelay = Instrument.instance().getRxTxDelay();
        gradientDelay = Instrument.instance().getGradDelay();
    }

    public double getGradientRxTxDelay() {
        return (rxTxDelay - gradientDelay);
    }
    public double getGradientDelay() {
        return (rxTxDelay - gradientDelay);
    }
    public double getRxTxDelay() {
        return rxTxDelay;
    }

}


