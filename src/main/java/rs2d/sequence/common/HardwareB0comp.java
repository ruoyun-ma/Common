package rs2d.sequence.common;

import rs2d.commons.log.Log;
import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.hardware.controller.gradient.GradientHandler;
import rs2d.spinlab.hardware.controller.peripherique.GradientHandlerInterface;
import rs2d.spinlab.hardware.devices.DeviceManager;

public class HardwareB0comp {
    private double amp = Double.NaN;
    private double phase = Double.NaN;
    private boolean status = false;

    public HardwareB0comp() {
        amp = 0;
        phase = 0;
        status = false;
        DeviceManager.getInstance().getGradientHandler().ifPresent(gradientHandler -> {
            if (gradientHandler.isAvailable() && gradientHandler.isConnected()) {
                try {
                    // available after V1.214
                    amp = gradientHandler.getB0compAmp().getValue().doubleValue();
                    phase = gradientHandler.getB0compPhase().getValue().doubleValue();
                    status = gradientHandler.isB0compEnabled();

                } catch (Exception e) {
                    System.out.println(" B0-Compensation not supported by this " );
                    Log.info(getClass(), " B0-Compensation not supported by this ");
                }
            }
        });
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