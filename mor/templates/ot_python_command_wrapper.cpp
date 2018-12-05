#include "Wrapper.h"
#include "OT.hxx"

#define WRAPPERNAME @FUNCTION_NAME@

// FUNC_INIT( WRAPPERNAME, {return WRAPPER_OK;} )

FUNC_INFO( WRAPPERNAME ,
           {
               p_info->inSize_ = getNumberOfVariables(p_exchangedData, WRAPPER_IN);
               p_info->outSize_ = getNumberOfVariables(p_exchangedData, WRAPPER_OUT);
               return WRAPPER_OK;
           } )

FUNC_EXEC( WRAPPERNAME, {
        struct WrapperExchangedData * p_exchangedData = CAST( struct WrapperExchangedData * , p_state ) ;
        int rc  = 0 ;
        char * temporaryDirectory = 0 ;
        char* cmd = 0 ;
        char * currentWorkingDirectory = 0 ;
        void * p_error = 0;

        currentWorkingDirectory = getCurrentWorkingDirectory(p_error);
        printf("CURRENT DIRECTORY = %s", currentWorkingDirectory);

        if (currentWorkingDirectory == NULL)
            return WRAPPER_EXECUTION_ERROR;

        if ( createInputFiles( currentWorkingDirectory, p_exchangedData , inPoint, p_error ) )
            return WRAPPER_EXECUTION_ERROR;

        rc = runInsulatedCommand( currentWorkingDirectory, p_exchangedData , inPoint, p_error );
        if( !rc )
            if ( readOutputFiles( currentWorkingDirectory, p_exchangedData , outPoint, p_error ) )
                return WRAPPER_EXECUTION_ERROR;

        free ( currentWorkingDirectory );
        return WRAPPER_OK;
    } )

FUNC_EXEC_SAMPLE( WRAPPERNAME, {

        int rc  = 0 ;
        char* cmd = 0 ;
        char * currentWorkingDirectory = 0 ;
        char * temporaryDirectory = 0 ;
        void * p_error = 0;

        currentWorkingDirectory = getCurrentWorkingDirectory(p_error);

        if (currentWorkingDirectory == NULL)
            return WRAPPER_EXECUTION_ERROR;

        const char FUNCTIONNAME[] = "FUNC_EXEC_SAMPLE_@FUNCTION_NAME@";
        printSample(FUNCTIONNAME, inSample);

        unsigned long size(inSample->size_);
        for (unsigned long i = 0; i < size; i++)
            {
                struct point mypointIn;
                struct point mypointOut;

                // Create inPoint
                mypointIn.size_ = inSample->dimension_;
                mypointIn.data_ = &( inSample->data_[i * inSample->dimension_] );
                const struct point mypointInConst = mypointIn;

                // Create outPoint
                mypointOut.size_ = outSample->dimension_;
                mypointOut.data_ = &( outSample->data_[i*outSample->dimension_] );

                if ( createInputFiles( currentWorkingDirectory, p_exchangedData , &mypointInConst, p_error ) )
                    return WRAPPER_EXECUTION_ERROR;

                rc = runInsulatedCommand( currentWorkingDirectory, p_exchangedData , &mypointInConst, p_error ) ;

                if( !rc )
                    if ( readOutputFiles( currentWorkingDirectory, p_exchangedData , &mypointOut, p_error ) )
                        return WRAPPER_EXECUTION_ERROR;
            }

        free ( currentWorkingDirectory );
        return WRAPPER_OK;
    } )
