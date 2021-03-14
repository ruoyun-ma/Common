package rs2d.sequence.common;

import rs2d.commons.log.Log;

public class ModelFactory {
    public ModelInterface getModel(String modelType, SeqPrep seqPrep){
        if(modelType == null){
            return null;
        }
        if(modelType.equalsIgnoreCase("ExtTrig")){
            System.out.println("xxx we create ExtTrig");
            return new ExtTrig(seqPrep);
        } else if(modelType.equalsIgnoreCase("FatSat")){
            System.out.println("xxx we create FatSat");
            return new FatSat(seqPrep);
        } else if(modelType.equalsIgnoreCase("FatSatWep")){
            System.out.println("xxx we create FatSatWep");
            return new FatSatWep(seqPrep);
        } else if(modelType.equalsIgnoreCase("SatBand")){
            System.out.println("xxx we create SatBand");
            return new SatBand(seqPrep);
        } else if(modelType.equalsIgnoreCase("TofSat")){
            System.out.println("xxx we create TofSat");
            return new TofSat(seqPrep);
        }
        Log.error(getClass(),"The required model is not supported!");

        return null;
    }
}
