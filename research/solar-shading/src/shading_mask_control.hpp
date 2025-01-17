



namespace Feel
{
    
    template <typename MeshType>
    bool 
    ShadingMask<MeshType>::testComparisonMaskValidationLevel1(
        std::string matrix_filename_NEW,std::string matrix_filename_CTRL)
    {
    
        bool QCTRL=false; 
        //std::string matrix_filename_NEW  = shadingMaskFolder+"/"+prefix_name+building_name+"_"+marker_name+".csv";
        //std::string matrix_filename_CTRL = shadingMaskFolder+"/"+prefix_name+"CTRL_"+building_name+"_"+marker_name+".csv";
        //cout<<"Test mat:"<<building_name+"_"+marker_name+".csv \n"; 
        std::ifstream FICH_NEW,FICH_CTRL;
        FICH_NEW.open(matrix_filename_NEW);
        FICH_CTRL.open(matrix_filename_CTRL);
        std::string MatLine_NEW,MatLine_CTRL;
        bool QRep=true; int i=0;
        while ((!FICH_CTRL.eof())) 
        {
            std::getline(FICH_NEW,MatLine_NEW);
            std::getline(FICH_CTRL,MatLine_CTRL);
            QRep=QRep && (MatLine_NEW.compare(MatLine_CTRL)==0); i++;
            if (!QRep) { QCTRL=true; }
        }
        FICH_NEW.close(); FICH_CTRL.close();
        //if (QCTRL) { cout<<"Test mat:"<<building_name+"_"+marker_name+".csv ===>"<<QRep<<"\n"; }
        return QRep;
    }


    template <typename MeshType>
    bool 
    ShadingMask<MeshType>::testComparisonMaskValidationLevel2(
        std::string matrix_filename_NEW,std::string matrix_filename_CTRL)
    {
        
        bool QCTRL=false;
        //std::string matrix_filename_NEW  = shadingMaskFolder+"/"+prefix_name+building_name+"_"+marker_name+".csv";
        //std::string matrix_filename_CTRL = shadingMaskFolder+"/"+prefix_name+"CTRL_"+building_name+"_"+marker_name+".csv";
        std::ifstream FICH_NEW,FICH_CTRL;
        FICH_NEW.open(matrix_filename_NEW);
        FICH_CTRL.open(matrix_filename_CTRL);
        std::string MatLine_NEW,MatLine_CTRL;
        bool QRep=true; int l=0;
        std::vector<std::vector<std::string> > parsedCsv_CTRL;
        std::vector<std::vector<std::string> > parsedCsv_NEW;

        int LiCTRL,LjCTRL;
        int LiNEW,LjNEW;

        double averageCTRL=0.0;
        double averageNEW=0.0;

        while ((!FICH_CTRL.eof())) 
        {   l++;
            std::getline(FICH_NEW,MatLine_NEW);
            std::getline(FICH_CTRL,MatLine_CTRL);
            std::stringstream lineStream_CTRL(MatLine_CTRL);
            std::stringstream lineStream_NEW(MatLine_NEW);
            std::string cell_CTRL;
            std::string cell_NEW;
            std::vector<std::string> parsedRow_CTRL;
            std::vector<std::string> parsedRow_NEW;
            while(std::getline(lineStream_CTRL,cell_CTRL,',')) { parsedRow_CTRL.push_back(cell_CTRL);}
            parsedCsv_CTRL.push_back(parsedRow_CTRL);

            while(std::getline(lineStream_NEW,cell_NEW,',')) { parsedRow_NEW.push_back(cell_NEW); }
            parsedCsv_NEW.push_back(parsedRow_NEW);
        }

        LjCTRL=parsedCsv_CTRL.size()-1; LiCTRL=parsedCsv_CTRL[0].size();
        LjNEW=parsedCsv_NEW.size()-1;   LiNEW=parsedCsv_NEW[0].size();
        for (int i=0; i<LiCTRL;i++) 
        {
                for (int j=0; j<LjCTRL;j++) {  averageCTRL=averageCTRL+strtod(parsedCsv_CTRL[j][i].c_str(), NULL);
                    //cout<<parsedCsv_CTRL[j][i]<<" "; 
                }
                //cout<<"\n";
                for (int j=0; j<LjCTRL;j++) {  averageNEW=averageNEW+strtod(parsedCsv_NEW[j][i].c_str(), NULL);
                    //cout<<parsedCsv_NEW[j][i]<<" "; 
                }
                //cout<<"\n";
        }

        averageCTRL=averageCTRL/(double(LiCTRL)*double(LjCTRL));
        averageNEW=averageNEW/(double(LiNEW)*double(LjNEW));

         
        if (QCTRL) {
            //cout<<"*** ERROR !!! ***\n";
            cout<<"  Name NEW  = ["<<matrix_filename_NEW<<"]\n";
            cout<<"  Name CTRL = ["<<matrix_filename_CTRL<<"]\n";    
            cout<<"  In matrix line "<<l<<"\n";  
            cout<<"  Lgn Size CTRL ="<<LiCTRL<<":"<<LjCTRL<<"\n";
            cout<<"  Lgn Size NEW  ="<<LiNEW<<":"<<LjNEW<<"\n"; 
            cout<<"  averageCTRL   ="<<std::fixed<<std::setprecision(9)<<averageCTRL<<"\n"; 
            cout<<"  averageNEW    ="<<std::fixed<<std::setprecision(9)<<averageNEW<<"\n"; 
        }
        FICH_NEW.close(); FICH_CTRL.close();
        //QCTRL=true;
        //if (QCTRL) { cout<<"Test mat:"<<building_name+"_"+marker_name+".csv ===>"<<QRep<<"\n"; 
        //            //getchar();
        //}
        return QRep;
        
       //return true;
    }


   template <typename MeshType>
    void 
    ShadingMask<MeshType>::testComparisonAllMasksValidation()
    {
        
        std::string prefix_name  = "SM_Matrix_";

        bool QTestOKLevel1=true;
        bool QTestOKLevel2=true;
        tic();
        cout<<"\n[INFO] : Control and validation  ===> ";
        if ((j_["/Buildings"_json_pointer].contains("aggregatedMarkers") ) && (true)) 
        {
            bool QTestOKLevel1=true;
            bool QTestOKLevel2=true;
            for(int i=0; i< M_listFaceMarkers.size(); i++)
            {
                std::string building_name = M_listFaceMarkers[i];
                std::string marker = "";
                std::string shadingMaskFolder = (boost::filesystem::path(Environment::appRepository())/("shadingMasks")).string();

                std::string matrix_filename_NEW  = shadingMaskFolder+"/"+prefix_name+building_name+"_"+marker+".csv";
                std::string matrix_filename_CTRL = shadingMaskFolder+"/"+prefix_name+"CTRL_"+building_name+"_"+marker+".csv";

                QTestOKLevel1=QTestOKLevel1 && (testComparisonMaskValidationLevel1(matrix_filename_NEW,matrix_filename_CTRL));
                QTestOKLevel2=QTestOKLevel2 && (testComparisonMaskValidationLevel2(matrix_filename_CTRL,matrix_filename_CTRL));
            }
            
        }
        if (QTestOKLevel1) { cout<<" >TRUE\n"; } else { cout<<" >ERROR\n"; }
        M_metadataJson["shadingMask"]["validation"]["Time"] = toc("Time CTRL Masks");
        M_metadataJson["shadingMask"]["validation"]["CTRL Level1"] = QTestOKLevel1;
        M_metadataJson["shadingMask"]["validation"]["CTRL Leval2"] = QTestOKLevel2;
        
    }   
}

