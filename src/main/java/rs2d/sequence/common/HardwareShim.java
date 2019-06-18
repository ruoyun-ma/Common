package rs2d.sequence.common;

import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.hardware.controller.peripherique.ShimHandlerInterface;
import rs2d.spinlab.hardware.controller.shim.ShimHandler;
import rs2d.spinlab.hardware.devices.DeviceManager;

import java.util.ArrayList;

//import rs2d.sequence.common.HardwareShim;
//         HardwareShim hardwareShim = new HardwareShim();
//         setParamValue(HARDWARE_SHIM, hardwareShim.getValue());
//         setParamValue(HARDWARE_SHIM_LABEL, hardwareShim.getLabel());

public class HardwareShim {
    private ArrayList<Number> shim = new ArrayList<>();
    private String shimText = new String();

    public HardwareShim() throws Exception {
        if (DeviceManager.getInstance().getShimHandler().isPresent()) {
            ShimHandler getShimHandler = DeviceManager.getInstance().getShimHandler().get();
            if (getShimHandler.isAvailable() && getShimHandler.isConnected()) {
                for (String param : getShimHandler.getAll()) {
                    if (shimText.length() == 0) {
                        shimText = new String(param);
                    } else {
                        shimText = shimText.concat(" ");
                        shimText = shimText.concat(param);
                    }
                    shim.add(getShimHandler.read(param).getValue());
                }
            } else {
                shimText = new String("NotConnected");
                shim.add(0);
            }
        }
    }

    public ArrayList<Number> getValue() {
        return shim;
    }

    public String getLabel() {
        return shimText;
    }
}


