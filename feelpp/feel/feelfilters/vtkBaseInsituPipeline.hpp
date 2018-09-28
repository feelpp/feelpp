
#ifndef FEELPP_VTKBASEINSITUPIPELINE_HPP
#define FEELPP_VTKBASEINSITUPIPELINE_HPP 1

#include <map>
#include <string>

#include <feel/feelcore/environment.hpp>

#if defined( __GNUC__ ) && !( defined( __clang__ ) )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkStdString.h>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPPipeline.h>
#include <vtkLiveInsituLink.h>
#include <vtkPVTrivialProducer.h>
#include <vtkSMProxy.h>
#include <vtkSMProxyManager.h>
#include <vtkSMSessionProxyManager.h>
#include <vtkSMSourceProxy.h>

namespace Feel
{

class vtkBaseInsituPipeline : public vtkCPPipeline
{
public:
    static vtkBaseInsituPipeline* New();
    vtkTypeMacro( vtkBaseInsituPipeline, vtkCPPipeline );
    virtual void PrintSelf( ostream& os, vtkIndent indent ) override;
    virtual void Initialize();

    virtual int RequestDataDescription( vtkCPDataDescription* dataDescription ) override;
    virtual int CoProcess( vtkCPDataDescription* dataDescription ) override;

protected:
    vtkBaseInsituPipeline();
    virtual ~vtkBaseInsituPipeline();

private:
    void CreatePipeline( vtkCPDataDescription* dataDescription, std::string inname );
    void UpdateProducers( vtkCPDataDescription* dataDescription );
    void DoLiveVisualization( vtkCPDataDescription* dataDescription );

    vtkBaseInsituPipeline( const vtkBaseInsituPipeline& ); // Not implemented
    void operator=( const vtkBaseInsituPipeline& );        // Not implemented
    int OutputFrequency;
    std::string FileName;

    bool pipelineCreated;

    vtkSMProxyManager* M_pxm;
    vtkSMSessionProxyManager* M_spxm;

    std::map<std::string, vtkSmartPointer<vtkSMSourceProxy>> M_producerMap;
    //std::map<std::string, vtkSmartPointer<vtkSMSourceProxy>> M_producerMap;

    vtkSmartPointer<vtkLiveInsituLink> M_insituLink;
    vtkSmartPointer<vtkSMProxy> M_proxy;
};
}

#if defined( __GNUC__ ) && !( defined( __clang__ ) )
#pragma GCC diagnostic pop
#endif

#endif
