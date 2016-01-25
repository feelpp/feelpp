
try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------
def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      dataProd = coprocessor.CreateProducer( datadescription, "input" )
      SetActiveSource(dataProd)

      # Prints data description
      # print str(datadescription.GetInputDescriptionByName("input").GetGrid())

      # create a new 'Parallel PolyData Writer'
      #parallelUnstructuredGridWriter0 = servermanager.writers.XMLPUnstructuredGridWriter(Input=dataProd)

      # register the writer with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the data, etc.
      #coprocessor.RegisterWriter(parallelUnstructuredGridWriter0, filename='filename_%t.pvtu', freq=1)

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  #freqs = {'input': [10, 100]}
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(True, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    hostname = "localhost"
    port = 22222

    userdata = datadescription.GetUserData()
    if(userdata != None):
        if( userdata.HasArray("hostname") ):
            hostname = userdata.GetAbstractArray("hostname").GetValue(0)
        if( userdata.HasArray("port") ):
            port = int(userdata.GetArray("port").GetTuple1(0))

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    #coprocessor.WriteData(datadescription);

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, hostname, port)
