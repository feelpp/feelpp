
#include <feel/feelfilters/vtkBaseInsituPipeline.hpp>

namespace Feel
{

#if defined(__GNUC__) && !(defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

vtkStandardNewMacro(vtkBaseInsituPipeline);

#if defined(__GNUC__) && !(defined(__clang__))
#pragma GCC diagnostic pop
#endif

//----------------------------------------------------------------------------
vtkBaseInsituPipeline::vtkBaseInsituPipeline()
{
    this->OutputFrequency = 0;
    this->pipelineCreated = false;

    /* Get current proxy manager */
    M_pxm = vtkSMProxyManager::GetProxyManager();
    M_spxm = M_pxm->GetActiveSessionProxyManager();
}
//----------------------------------------------------------------------------
vtkBaseInsituPipeline::~vtkBaseInsituPipeline()
{
}
//----------------------------------------------------------------------------
void vtkBaseInsituPipeline::Initialize(/*int outputFrequency, std::string& fileName*/)
{
    //this->OutputFrequency = outputFrequency;
    //this->FileName = fileName;
}
//----------------------------------------------------------------------------
int vtkBaseInsituPipeline::RequestDataDescription(
vtkCPDataDescription* dataDescription)
{
    if(!dataDescription)
    {
        vtkWarningMacro("dataDescription is NULL.");
        return 0;
    }

    /*
    if(dataDescription->GetForceOutput() == true)
    {
        for(int i = 0; i < dataDescription->GetNumberOfInputDescriptions(); i++)
        {
            dataDescription->GetInputDescription(i)->AllFieldsOn();
            dataDescription->GetInputDescription(i)->GenerateMeshOn();
        }
        return 1;
    }
    */

    int timestep = dataDescription->GetTimeStep();
    int ninputs = dataDescription->GetNumberOfInputDescriptions();

    for(int i = 0; i < ninputs; i++)
    {
        dataDescription->GetInputDescription(i)->AllFieldsOn();
        dataDescription->GetInputDescription(i)->GenerateMeshOn();
    }

    return 1;
}

void vtkBaseInsituPipeline::CreatePipeline(vtkCPDataDescription* dataDescription, std::string inname)
{
    if(!(dataDescription->GetInputDescriptionByName(inname.c_str())))
    {
        std::cout << "Simulation input name " << inname.c_str() << "does not exists" << std::endl;
        return;
    } 

    /* Get the data to output on ParaView */
    vtkSmartPointer<vtkMultiBlockDataSet> grid = vtkMultiBlockDataSet::SafeDownCast(
            dataDescription->GetInputDescriptionByName(inname.c_str())->GetGrid());

    // Create a vtkPVTrivialProducer and set its output
    // to be the input grid.
    vtkSmartPointer<vtkSMSourceProxy> producer;

    /* Create a proxy for the data producer that will be created in paraview interface */
    vtkSMProxy * px = M_spxm->NewProxy("sources", "PVTrivialProducer");
    producer.TakeReference(vtkSMSourceProxy::SafeDownCast(px));
    /* It is very important to register the proxy */
    /* The NewProxy previously called doesn't do it */
    /* Ref: Qt/Core/pqObjectBuilder.cxx */
    M_spxm->RegisterProxy("sources", producer);

    /* Get the client side object: The object in ParaView */
    /* and set its output */
    vtkObjectBase* clientSideObject = producer->GetClientSideObject();
    vtkPVTrivialProducer* realProducer = vtkPVTrivialProducer::SafeDownCast(clientSideObject);
    realProducer->SetOutput(grid);

    /* Update the producer */
    producer->UpdatePipeline();

    /* Record the producer for updates */
    M_producerMap[inname.c_str()] = producer;
}

void vtkBaseInsituPipeline::UpdateProducers(vtkCPDataDescription* dataDescription)
{
    /* Ensure that the producers are created */
    if(!pipelineCreated)
    {
        this->CreatePipeline(dataDescription, std::string("input"));
        this->pipelineCreated = true;
    }

    /* Update the output of each recorded producer */
    for(auto it = this->M_producerMap.begin(); it != this->M_producerMap.end(); it++)
    {
        vtkObjectBase* clientSideObject = it->second->GetClientSideObject();
        vtkPVTrivialProducer* realProducer = vtkPVTrivialProducer::SafeDownCast(clientSideObject);
        realProducer->SetOutput(dataDescription->GetInputDescriptionByName(it->first.c_str())->GetGrid(), dataDescription->GetTime());
    }
}

void vtkBaseInsituPipeline::DoLiveVisualization(vtkCPDataDescription* dataDescription)
{
    /* Initialize the insitu link, if not already done */
    if(!M_insituLink)
    {
        M_insituLink = vtkSmartPointer<vtkLiveInsituLink>::New();
        M_insituLink->SetHostname( soption( _name="exporter.vtk.insitu.hostname" ).c_str() );
        M_insituLink->SetInsituPort( ioption( _name="exporter.vtk.insitu.port" ) );
        M_insituLink->Initialize( M_spxm );  
    }

    /* Get the time information from the data */
    double time = dataDescription->GetTime();
    vtkIdType timeStep = dataDescription->GetTimeStep();

    /* Go into the main loop for visualization */
    while(true)
    {
        /* Initialize the link and update data */
        M_insituLink->InsituUpdate(time, timeStep); 
        for(auto it = M_producerMap.begin(); it != M_producerMap.end(); it++)
        {
            /* Update the pipeline for the current timestep */ 
            it->second->UpdatePipeline(time);
        }
        M_insituLink->InsituPostProcess(time, timeStep); 

        /* Wait for changes on ParaView side (Pause ...) */
        if(M_insituLink->GetSimulationPaused())
        {
            if(M_insituLink->WaitForLiveChange())
            {
                break;
            }
        }
        else
        {
            break;
        }
    }
}

//----------------------------------------------------------------------------
int vtkBaseInsituPipeline::CoProcess(
vtkCPDataDescription* dataDescription)
{
    if(!dataDescription)
    {
        vtkWarningMacro("DataDescription is NULL");
        return 0;
    }

    /* Get the grid to output in ParaView */
    vtkSmartPointer<vtkMultiBlockDataSet> grid = vtkMultiBlockDataSet::SafeDownCast(
            dataDescription->GetInputDescriptionByName("input")->GetGrid());
    if(grid == NULL)
    {
        vtkWarningMacro("DataDescription is missing input unstructured grid.");
        return 0;
    }
    if(this->RequestDataDescription(dataDescription) == 0)
    {
        return 1;
    }

    /* see coprocessing.py */

    /* Update the producers and do live visualization */
    this->UpdateProducers(dataDescription);
    this->DoLiveVisualization(dataDescription);

    return 1;
}
//----------------------------------------------------------------------------
void vtkBaseInsituPipeline::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    /*
    os << indent << "OutputFrequency: " << this->OutputFrequency << "\n";
    os << indent << "FileName: " << this->FileName << "\n";
    */
}

} // namespace Feel
