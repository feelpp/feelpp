#ifndef FEEL_OMMODEL_HPP
#define FEEL_OMMODEL_HPP

namespace Feel
{


class OMModel
{
public :
    OMModel()
    {}

    int run( int argc, char** argv )
    {
        int res;
        DATA data;
        MODEL_DATA modelData;
        SIMULATION_INFO simInfo;
        data.modelData = &modelData;
        data.simulationInfo = &simInfo;
        measure_time_flag = 0;
        compiledInDAEMode = 0;
        MMC_INIT(0);
        omc_alloc_interface.init();
        {
            MMC_TRY_TOP()
                MMC_TRY_STACK()

                this->setupDataStruc(&data, threadData);
            res = _main_SimulationRuntime(argc, argv, &data, threadData);

            MMC_ELSE()
                rml_execution_failed();
            fprintf(stderr, "Stack overflow detected and was not caught.\nSend us a bug report at https://trac.openmodelica.org/OpenModelica/newticket\n    Include the following trace:\n");
            printStacktraceMessages();
            fflush(NULL);
            return 1;
            MMC_CATCH_STACK()
                MMC_CATCH_TOP(return rml_execution_failed());
        }

        fflush(NULL);
        EXIT(res);
        return res;
    }

    virtual void setupDataStruc( DATA *data, threadData_t *threadData )=0;


private :

};


}//namespace Feel
#endif
