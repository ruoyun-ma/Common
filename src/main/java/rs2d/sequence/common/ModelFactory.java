package rs2d.sequence.common;

import rs2d.commons.log.Log;

/**
 *  V1.0 - 2021.3 XG
 *
 */

public class ModelFactory {
    public ModelInterface getModel(String modelType, SeqPrep seqPrep){
        if(modelType == null){
            return null;
        }
        if(modelType.equalsIgnoreCase("ExtTrig")){
            return new ExtTrig(seqPrep);
        } else if(modelType.equalsIgnoreCase("FatSat")){
            return new FatSat(seqPrep);
        } else if(modelType.equalsIgnoreCase("FatSatWep")){
            return new FatSatWep(seqPrep);
        } else if(modelType.equalsIgnoreCase("SatBand")){
            return new SatBand(seqPrep);
        } else if(modelType.equalsIgnoreCase("TofSat")){
            return new TofSat(seqPrep);
        } else if(modelType.equalsIgnoreCase("FlowComp")){
            return new FlowComp(seqPrep);
        } else if(modelType.equalsIgnoreCase("InvRec")){
            return new InvRec(seqPrep);
        } else if(modelType.equalsIgnoreCase("VFL")){
            return new VFL(seqPrep);
        }

        Log.error(getClass(),"The required model is not supported!");

        return null;
    }
}
