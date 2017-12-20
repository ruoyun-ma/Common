package rs2d.sequence.common;

import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.hardware.controller.peripherique.GradientHandlerInterface;

public class HardwareB0comp {
    private double amp = Double.NaN;
    private double phase = Double.NaN;
    private boolean status = false;

    public HardwareB0comp() {
        GradientHandlerInterface gradientHandler = HardwareHandler.getInstance().getGradientHandler();
        amp = 0;
        phase = 0;
        status = false;
        if (gradientHandler.isAvailable() && gradientHandler.isConnected()) {
            try {
                amp = gradientHandler.getB0compAmp().getValue().doubleValue();
                phase = gradientHandler.getB0compPhase().getValue().doubleValue();
                status = gradientHandler.isB0compEnabled();
            } catch (Exception e) {
                System.out.println(" B0-Compensation not supported by this ");
            }
        }
    }

    public double getAmp() {
        return amp;
    }

    public boolean getStatus() {
        return status;
    }

    public double getB0compPhase() {
        return phase;
    }
}