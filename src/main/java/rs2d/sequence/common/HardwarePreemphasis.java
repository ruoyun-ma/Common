package rs2d.sequence.common;

import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.hardware.controller.peripherique.GradientHandlerInterface;

import java.util.ArrayList;

public class HardwarePreemphasis {
    private ArrayList<Number> preemphasisT = new ArrayList<>();
    private ArrayList<Number> preemphasisA = new ArrayList<>();
    private ArrayList<Number> A0 = new ArrayList<>();
    private ArrayList<Number> DC = new ArrayList<>();

    public HardwarePreemphasis() throws Exception {
        while (preemphasisT.size() < 18) {
            preemphasisT.add(0);
        }
        while (preemphasisA.size() < 18) {
            preemphasisA.add(0);
        }
        while (DC.size() < 4) {
            DC.add(0);
        }
        while (A0.size() < 3) {
            A0.add(0);
        }

        GradientHandlerInterface gradientHandler = HardwareHandler.getInstance().getGradientHandler();
        if (gradientHandler.isAvailable() && gradientHandler.isConnected()) {
            for (String param : gradientHandler.getAll()) {
                int v = (param.charAt(0) == 'x') ? 0 : ((param.charAt(0) == 'y') ? 1 : ((param.charAt(0) == 'z') ? 2 : 3));
                if (param.charAt(1) == 'T' && param.charAt(2) != '0') {                 // T1 T2 T3
                    v = v * 3 + Integer.parseInt(String.valueOf(param.charAt(2))) - 1;
                    preemphasisT.set(v, gradientHandler.read(param).getValue());

                } else if (param.charAt(1) == 'A' && param.charAt(2) != '0') {          // A1 A2 A3
                    v = v * 3 + Integer.parseInt(String.valueOf(param.charAt(2))) - 1;
                    preemphasisA.set(v, gradientHandler.read(param).getValue());

                } else if (param.charAt(1) == 'A' && param.charAt(2) == '0') {          // A0
                    A0.set(v, gradientHandler.read(param).getValue());

                } else if (param.charAt(1) == 'D') {                                    // DC
                    DC.set(v, gradientHandler.read(param).getValue());

                } else if (param.charAt(0) == 'B') {                                    // DC
                    DC.set(v, gradientHandler.read(param).getValue());

                } else if (param.charAt(1) == 'B' && param.charAt(3) == 'T') {          // B0T1 B0T2 B0T3
                    v = 9 + v * 3 + Integer.parseInt(String.valueOf(param.charAt(4))) - 1;
                    preemphasisT.set(v, gradientHandler.read(param).getValue());

                } else if (param.charAt(1) == 'B' && param.charAt(3) == 'A') {           // B0A1 B0A2 B0A3
                    v = 9 + v * 3 + Integer.parseInt(String.valueOf(param.charAt(4))) - 1;
                    preemphasisA.set(v, gradientHandler.read(param).getValue());
                }
            }
        }
    }

    public ArrayList<Number> getA0() {
        return A0;
    }

    public ArrayList<Number> getDC() {
        return DC;
    }

    public ArrayList<Number> getTime() {
        return preemphasisT;
    }

    public ArrayList<Number> getAmplitude() {
        return preemphasisA;
    }

}


