package rs2d.sequence.common;

import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.hardware.controller.peripherique.ShimHandlerInterface;

import java.util.ArrayList;

public class HardwareShim {
    private ArrayList<Number> shim = new ArrayList<>();
    private String shimText = new String();

    public HardwareShim() throws Exception {
        ShimHandlerInterface getShimHandler = HardwareHandler.getInstance().getShimHandler();
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

    public ArrayList<Number> getValue() {
        return shim;
    }

    public String getLabel() {
        return shimText;
    }
}


