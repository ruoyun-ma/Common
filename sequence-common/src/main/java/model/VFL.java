package model;

import common.RFPulse;
import kernel.SeqPrep;
import org.apache.commons.math3.complex.Complex;
import rs2d.commons.log.Log;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.table.Order;

import java.util.ArrayList;
import java.util.Collections;

import static common.CommonUP.*;

/**
 * V1.0 - 2021.3 XG
 */

public class VFL implements ModelInterface {
    public SeqPrep parent;
    protected static boolean isVFLEnabled;

    public int echoTrainLengthEncoded;
    public boolean isLongVFL;
    private PSS pss = new PSS();
    private EPG epg = new EPG();
    private Order LoopOrder = Order.LoopB;

    private double T1 = 1.0;
    private double T2 = 0.1;
    private ArrayList<Double> sigCoherence; // without relaxation
    private ArrayList<Double> sigEstimate;
    private ArrayList<Double> variableFA;

    protected enum UP implements GeneratorParamEnum {
        FLIP_ANGLE_REFOC_VFL,
        ECHO_TRAIN_LENGTH_ENCODED,
        VFL_FLIP_ANGLE_TRAIN,
        FLIP_ANGLE_MIN_VFL,
        FLIP_ANGLE_MID_VFL,
        ;

        @Override
        public Param build() {
            throw new UnsupportedOperationException("Unable to create a Param");
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
    }

    public VFL() {
    }

    public VFL(SeqPrep parent) {
        this.parent = parent;
    }

    public void init() {
        isVFLEnabled = parent.getBoolean(UP.FLIP_ANGLE_REFOC_VFL);
    }

    @Override
    public void init(SeqPrep parent) {
        this.parent = parent;
        init();
    }

    @Override
    public void initFinal() throws Exception {
        iniETLandVFL();
    }

    @Override
    public void prep() throws Exception {
        clcPulse(); //just in case the user skip this step in the sequence
    }

    @Override
    public void prepFinal() {
        Log.warning(getClass(), ("The order of vfl is LoopOrder = " + LoopOrder.name()));
        parent.pulseTX180.prepTxAmpMultiFA(parent.getListInt(TX_ROUTE), parent.getListDouble(UP.VFL_FLIP_ANGLE_TRAIN), LoopOrder);
    }

    @Override
    public RFPulse getRfPulses() {
        if (parent.pulseTX180 != null) {
            return parent.pulseTX180;
        } else {
            return null;
        }
    }

    @Override
    public double getDuration() {
        return 0;
    }

    @Override
    public String getName() {
        return "VFL";
    }

    @Override
    public boolean isEnabled() {
        return isVFLEnabled;
    }

    protected void clcPulse() {
        if (!parent.pulseTX180.checkPower(180, parent.observeFrequency, parent.nucleus)) {
            parent.notifyOutOfRangeParam(TX_LENGTH_180, parent.pulseTX180.getPulseDuration(), ((NumberParam) parent.getParam(TX_LENGTH_180)).getMaxValue(), "Pulse length too short to reach RF power with this pulse shape");
            parent.getParam(TX_LENGTH_180).setValue(parent.pulseTX180.getPulseDuration());
        }
    }

    public ArrayList<Double> getVFL() {
        if (isVFLEnabled) {
            if (isLongVFL) {
                return variableFA;
            } else {
                ArrayList<Double> pseudoVFA = new ArrayList<>(Collections.nCopies(parent.getInt(ECHO_TRAIN_LENGTH), parent.getDouble(FLIP_ANGLE_REFOC)));
                pseudoVFA.set(0, parent.getDouble(FLIP_ANGLE) + parent.getDouble(FLIP_ANGLE_REFOC) / 2);
                return pseudoVFA;
            }
        } else {
            return new ArrayList<>(Collections.nCopies(parent.getInt(ECHO_TRAIN_LENGTH), parent.getDouble(FLIP_ANGLE_REFOC)));
        }
    }

    public int getNumPSSPrep() {
        return pss.numPSSPrep;
    }

    public double getEffTE() {
        if (isLongVFL)
            return getEffTE(parent.getDouble(ECHO_TIME), parent.getDouble(ECHO_SPACING), parent.getInt(ECHO_TRAIN_LENGTH));
        else
            return parent.getDouble(ECHO_TIME);
    }

    public double getInvEffTE() {
        return getInvEffTE(parent.getDouble(ECHO_TIME), parent.getDouble(ECHO_SPACING), parent.getInt(ECHO_TRAIN_LENGTH));
    }

    public void setOrder(Order order) {
        // we provide an API for user input order
        LoopOrder = order;
    }

    private void iniETLandVFL() throws Exception {
        double fa_min = parent.getDouble(UP.FLIP_ANGLE_MIN_VFL);
        double fa_mid = parent.getDouble(UP.FLIP_ANGLE_MID_VFL);

        if (fa_min <= 0.0)
            fa_min = Math.max(parent.getDouble(FLIP_ANGLE_REFOC) / 6.0 + 5.0, 20.0);
        if (fa_mid <= 0.0)
            fa_mid = parent.getDouble(FLIP_ANGLE_REFOC) / 2.0 + 10.0;

        if (isVFLEnabled && parent.echoTrainLength >= 20 && parent.getDouble(FLIP_ANGLE_REFOC) < 180.0) {
            echoTrainLengthEncoded = 1;
            while (Math.floorMod(echoTrainLengthEncoded, 2) != 0) { //we need a even number of echoTrainLengthEncoded //Maybe NOT...
                getVFL(fa_min, fa_mid, parent.getDouble(FLIP_ANGLE_REFOC), parent.echoTrainLength); // our method is independent of ESP
                echoTrainLengthEncoded = parent.echoTrainLength - getNumPSSPrep();  //encoded echoTrainLength shrinks
                System.out.println("echoTrainLength " + parent.echoTrainLength + " echoTrainLengthScan " + echoTrainLengthEncoded + " + numPSSPrep " + getNumPSSPrep());
                parent.echoTrainLength++;

                if (echoTrainLengthEncoded < 1) {
                    Log.error(getClass(), "Encoding Matrix ERROR : number of echoTrainLength too short for variable flip angle scheme\n");
                }
            }
            parent.echoTrainLength = echoTrainLengthEncoded + getNumPSSPrep();
        } else {
            echoTrainLengthEncoded = parent.echoTrainLength;
        }

        parent.getParam(UP.VFL_FLIP_ANGLE_TRAIN).setValue(getVFL());
        parent.getParam(ECHO_TRAIN_LENGTH).setValue(parent.echoTrainLength);
        parent.getParam(UP.ECHO_TRAIN_LENGTH_ENCODED).setValue(echoTrainLengthEncoded);
    }

    private double[] getVFL(double flipAngleMin, double flipAnlgeCenter, double flipAngleMax, int ETL) {
        // Reed Busse, MRM, 2008
        ArrayList<Double> sigTargetPSS;
        this.variableFA = new ArrayList<>();

        double sigFirst = Math.pow(Math.sin(Math.toRadians(flipAngleMax) / 2.0), 2);
        double sigMin = this.pss.getPSSControl(flipAngleMin);
        double sigCenter = this.pss.getPSSControl(flipAnlgeCenter);
        double sigMax = this.pss.getPSSControl(flipAngleMax);

        sigTargetPSS = this.pss.getPSS(sigFirst, sigMin, sigCenter, sigMax, ETL, 0.99);

        this.sigCoherence = this.epg.cpmg_epg(ETL, 90, 0.0, 0.0, 1.0, flipAngleMax, sigTargetPSS, this.variableFA); //ESP is a fake input
        //this.sigEstimate = this.epg.cpmg_epg(ETL, 90, T1, T2, ESP, flipAngleMax, sigTargetPSS, this.variableFA);

        //System.out.println("!!! PSS Target: " + sigTargetPSS.toString());
        //System.out.println("!!! PSS Coherence: " + sigCoherence.toString());
        //System.out.println("!!! PSS Estimate: " + sigEstimate.toString());

        if (!this.variableFA.isEmpty()) {
            isLongVFL = true;
        } else {
            isLongVFL = false;
        }

        return variableFA.stream().mapToDouble(Double::doubleValue).toArray();
        //return sigEstimate.stream().mapToDouble(Double::doubleValue).toArray();
    }

    private double getEffTE(double physical_te, double ESP, int ETL) {
        ArrayList<Double> arrayList = new ArrayList<>();
        this.sigEstimate = this.epg.cpmg_epg(ETL, 90, T1, T2, ESP, .0, new ArrayList<>(), this.variableFA);

        if (this.sigEstimate != null && !this.sigEstimate.isEmpty() && this.sigCoherence != null && !this.sigCoherence.isEmpty()) {
            for (int i = 0; i < this.sigEstimate.size(); i++) {
                arrayList.add(i, -this.T2 * Math.log(this.sigEstimate.get(i) / this.sigCoherence.get(i)));
            }
        }
        //System.out.println("getEffTE "+arrayList.toString());

        return arrayList.get((int) Math.round(physical_te / ESP) - 1);
    }

    private double getInvEffTE(double contrast_te, double ESP, int ETL) {
        ArrayList<Double> arrayList = new ArrayList<>();
        this.sigEstimate = this.epg.cpmg_epg(ETL, 90, T1, T2, ESP, .0, new ArrayList<>(), this.variableFA);

        if (this.sigEstimate != null && !this.sigEstimate.isEmpty() && this.sigCoherence != null && !this.sigCoherence.isEmpty()) {
            for (int i = 0; i < this.sigEstimate.size(); i++) {
                arrayList.add(Math.abs(this.sigEstimate.get(i) - this.sigCoherence.get(i) * Math.exp(-contrast_te / this.T2)));
            }
        }
        return arrayList.indexOf(Collections.min(arrayList)) * ESP;
    }
}

class PSS {
    protected int numPSSPrep = 0;

    protected PSS() {
    }

    protected double getPSSControl(double flipAngle) {
        // David Alsop, MRM, 1997
        double Sig_PSS = 0.0;
        int Nspins = 2000;

        double term = 0.5 * Math.pow(Math.sin(Math.toRadians(flipAngle)) / (1 - Math.cos(Math.toRadians(flipAngle))), 2);
        //System.out.println("term: "+ term);
        for (int dPhi = 0; dPhi < Nspins + 1; dPhi++) {
            Sig_PSS += Math.pow(1 + term - term * Math.cos((dPhi / (double) Nspins) * 2.0 * Math.PI), -0.5);
        }
        System.out.println("PSS Signal Control Point: " + flipAngle + " | " + Sig_PSS / Nspins);
        return Sig_PSS /= Nspins;
    }

    protected ArrayList<Double> getPSS(double sigFirst, double sigMin, double sigCenter, double sigMax, int ETL, double PSSPrepQuality) {
        ArrayList<Double> sigTargetPSS = new ArrayList<>();

        double sigSeg1 = sigFirst;
        while ((sigSeg1 - sigMin) / sigMin > 1 - PSSPrepQuality) {
            sigTargetPSS.add(sigSeg1);
            sigSeg1 = 0.5 * (sigSeg1 + sigMin);
            //System.out.println("sigSeg1 "+sigSeg1+" sigMin "+sigMin+" diff = "+ (sigSeg1-sigMin)/sigMin);
        }

        // Three segments
        int seg1 = sigTargetPSS.size();
        this.numPSSPrep = seg1;

        int seg2 = (int) Math.floor(ETL / 2.0) - seg1;
        int seg3 = (int) Math.ceil(ETL / 2.0);

        for (int i = 0; i < seg2; i++) {
            sigTargetPSS.add(sigMin + (sigCenter - sigMin) * Math.pow((i + 1) / (double) seg2, 1.5));
        }
        for (int i = 0; i < seg3; i++) {
            sigTargetPSS.add(sigCenter + (sigMax - sigCenter) * Math.pow((i + 1) / (double) seg3, 1));
        }
        //System.out.println(sigTargetPSS);
        return sigTargetPSS;
    }
}

class EPG {
    protected EPG() {
    }

    protected ArrayList<Complex> rf_pulse_matrix(double alpha, double phi) {
        ArrayList<Complex> rfpulse = new ArrayList<>(9);

        double alpha_half = Math.toRadians(alpha) / 2.0;
        double c_a_h = Math.cos(alpha_half);
        double c_a_h2 = c_a_h * c_a_h;
        double s_a_h = Math.sin(alpha_half);
        double s_a_h2 = s_a_h * s_a_h;
        double c_a = Math.cos(Math.toRadians(alpha));
        double s_a = Math.sin(Math.toRadians(alpha));
        Complex i_phi = new Complex(0.0, phi);
        Complex i_half = new Complex(0.0, 0.5);

        rfpulse.add(0, new Complex(c_a_h2));
        rfpulse.add(1, new Complex(s_a_h2).multiply(i_phi.multiply(2.0).exp()));
        rfpulse.add(2, new Complex(s_a).multiply(i_phi.exp()).multiply(new Complex(0.0, -1.0)));

        rfpulse.add(3 + 0, new Complex(s_a_h2).multiply(i_phi.multiply(-2.0).exp()));
        rfpulse.add(3 + 1, new Complex(c_a_h2));
        rfpulse.add(3 + 2, new Complex(s_a).multiply(i_phi.multiply(-1.0).exp()).multiply(new Complex(0.0, 1.0)));

        rfpulse.add(6 + 0, new Complex(s_a).multiply(i_phi.multiply(-1.0).exp()).multiply(i_half).multiply(-1.0));
        rfpulse.add(6 + 1, new Complex(s_a).multiply(i_phi.multiply(1.0).exp()).multiply(i_half).multiply(1.0));
        rfpulse.add(6 + 2, new Complex(c_a));

        return rfpulse;
    }

    protected int mult_pulse_FpFmZ(ArrayList<Complex> rfpulse, ArrayList<Complex> FpFmZ, int num_cols_FpFmZ) {
        // This function does the matrix multiplication of rfpulse x FpFmZ
        // The result is returned in FpFmZ
        // rfpulse is a 2-D matrix of 3 rows and 3 columns represented as a 1-D array
        // FpFmZ is a 2-Dmatrix of 3 rows and num_cols_FpFmZ columns represented as a 1-D array

        final int num_rows = 3;
        final int num_cols_rfpulse = 3;

        ArrayList<Complex> result = new ArrayList<>(num_rows * num_cols_FpFmZ);
        Complex sum = new Complex(0.0);

        for (int r = 0; r < num_rows; r++) // row of rfpulse and result
        {
            for (int c = 0; c < num_cols_FpFmZ; c++) //columns of result and FpFmZ
            {
                for (int k = 0; k < num_rows; k++) // row index of FpFmZ
                {
                    sum = sum.add(rfpulse.get(r * num_cols_rfpulse + k).multiply(FpFmZ.get(k * num_cols_FpFmZ + c)));
                }
                result.add(r * num_cols_FpFmZ + c, sum);
                sum = new Complex(0.0);  // reset sum to zero
            }
        }
        for (int i = 0; i < num_rows * num_cols_FpFmZ; i++)
            FpFmZ.set(i, result.get(i));

        return 0;
    }

    protected int mult_ee_FpFmZ(ArrayList<Double> ee, ArrayList<Complex> FpFmZ, int num_cols_FpFmZ) {
        // This function does the matrix multiplication of ee x FpFmZ
        // The result is returned in FpFmZ
        // ee is a 2-D matrix of 3 rows and 3 columns represented as a 1-D array
        // FpFmZ is a 2-Dmatrix of 3 rows and num_cols_FpFmZ columns represented as a 1-D array

        final int num_rows = 3;
        final int num_cols_rfpulse = 3;

        ArrayList<Complex> result = new ArrayList<>(3 * num_cols_FpFmZ);
        Complex sum = new Complex(0.0);

        for (int r = 0; r < num_rows; r++) {
            for (int c = 0; c < num_cols_FpFmZ; c++) {
                for (int k = 0; k < num_rows; k++) {
                    sum = sum.add(FpFmZ.get(k * num_cols_FpFmZ + c).multiply(ee.get(r * num_cols_rfpulse + k)));
                }
                result.add(r * num_cols_FpFmZ + c, sum);
                sum = new Complex(0.0);
            }
        }
        for (int i = 0; i < num_rows * num_cols_FpFmZ; i++)
            FpFmZ.set(i, result.get(i));

        return 0;
    }

    protected void epg_grad(ArrayList<Complex> FpFmZ, int num_cols) {
        for (int i = 0; i < num_cols - 1; i++)
            FpFmZ.set(num_cols - 1 - i, FpFmZ.get(num_cols - 1 - i - 1));

        FpFmZ.set(0, new Complex(0.0));

        for (int i = 0; i < num_cols - 1; i++)
            FpFmZ.set(num_cols + i, FpFmZ.get(num_cols + i + 1));

        FpFmZ.set(num_cols + num_cols - 1, new Complex(0.0));
        FpFmZ.set(0, FpFmZ.get(num_cols).conjugate());
    }

    protected double getFlipAngleBasedonTargetSignal(ArrayList<Complex> FpFmZ, int ncols, double sigTarget, double T2, double ESP, double faMax) {
        Complex Z1 = FpFmZ.get(2 * ncols + 1);
        Complex Fp1 = FpFmZ.get(1);
        Complex Fm1 = FpFmZ.get(ncols + 1);

        double termP1, termM1;
        if (T2 == 0.0) {
            termP1 = Fp1.getReal() - sigTarget;
            termM1 = Fm1.getReal() - sigTarget;
        } else {
            termP1 = Fp1.getReal() - sigTarget * Math.exp(ESP / 2.0 / T2);
            termM1 = Fm1.getReal() - sigTarget * Math.exp(ESP / 2.0 / T2);
        }

        double term = Math.pow(Z1.abs(), 2) - termP1 * termM1;
        if (term < 0.0)
            return faMax;
        else if (termP1 == 0.0) {
            return 90.0;
        } else {
            return Math.abs(Math.toDegrees(2 * Math.atan((Z1.abs() - Math.pow(term, 0.5)) / termP1)));
        }
    }

    protected ArrayList<Double> cpmg_epg(int Nechos, double rf_90, double T1, double T2, double ESP, double faMax, ArrayList<Double> sigTarget, ArrayList<Double> variableFA) {
        // Nechos : number of echos to calculate
        // rf_90    : rf pulse value for initial excitation given in degrees
        // T1       : value for T1 relaxation time given in ms
        // T2       : value for T2 relaxation time given in ms
        // ESP      : value for CPMG echo spacing given in ms
        // sigTarget: value for the target PSS signal
        // return   :  signal : final calculated EPG signal, Nechos in length

        int nrows = 3;  // Fp, Fm, Z
        int ncols = 2 * Nechos;

        ArrayList<Complex> FpFmZ = new ArrayList<>(Collections.nCopies(nrows * ncols, new Complex(0.0)));
        ArrayList<Double> ee = new ArrayList<>(Collections.nCopies(9, 0.0));
        ArrayList<Complex> r90;
        ArrayList<Complex> r180;
        ArrayList<Double> signal = new ArrayList<>();

        // ****************************************
        // Initialize diagonal of relaxation matrix
        // ****************************************

        if (T2 == 0.0) {
            ee.set(0, 1.0);
            ee.set(4, 1.0);
        } else {
            ee.set(0, Math.exp(-ESP / 2.0 / T2));
            ee.set(4, Math.exp(-ESP / 2.0 / T2));
        }
        if (T1 == 0.0) {
            ee.set(8, 1.0);
        } else {
            ee.set(8, Math.exp(-ESP / 2.0 / T1));
        }
        //System.out.println("ee: "+ ee.toString());

        // ********************************
        // Set Z magnetization to 1.0 +0.0j
        // ********************************

        FpFmZ.set(2 * ncols + 0, new Complex(1.0));

        // *********************************
        // Create 90(y) rotation matrix
        // *********************************

        r90 = rf_pulse_matrix(rf_90, Math.PI / 2);

        // *********************************
        // Apply 90 degree pulse
        // *********************************

        mult_pulse_FpFmZ(r90, FpFmZ, ncols);

        // ****************************************
        // Apply 180 pulses of CPMG sequence
        // ****************************************

        for (int i = 0; i < Nechos; i++) {
            //Apply relaxation for  first half of half echo time
            mult_ee_FpFmZ(ee, FpFmZ, ncols);
            FpFmZ.set(2 * ncols, FpFmZ.get(2 * ncols).add(1.0 - ee.get(8)));

            //Apply gradient
            epg_grad(FpFmZ, ncols);

            // Create 180(x) rotation matrix
            double rf_180;
            if (T1 * T2 == 0.0) {
                rf_180 = getFlipAngleBasedonTargetSignal(FpFmZ, ncols, sigTarget.get(i), T2, ESP, faMax);
            } else {
                rf_180 = variableFA.get(i);
            }

            r180 = rf_pulse_matrix(rf_180, 0.0);

            // Apply effective 180 degree pulse
            mult_pulse_FpFmZ(r180, FpFmZ, ncols);

            //Apply relaxation for  second half of half echo time
            mult_ee_FpFmZ(ee, FpFmZ, ncols);
            FpFmZ.set(2 * ncols, FpFmZ.get(2 * ncols).add(1.0 - ee.get(8)));

            //Apply gradient
            epg_grad(FpFmZ, ncols);

            //System.out.println("Echo "+i+" FpFmZ: "+ FpFmZ.toString());
            //Save signal
            signal.add(i, FpFmZ.get(0).getReal());
            if (T1 * T2 == 0.0) {
                variableFA.add(rf_180);
            }
        }
        return signal;
    }

}