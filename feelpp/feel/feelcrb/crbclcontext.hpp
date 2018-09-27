/* -*- mode: c++ -*-
 */
/**
   \file crbCLContext.hpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2014-05-01
 */
#ifndef __CRBCLCONTEXT_H
#define __CRBCLCONTEXT_H 1

#ifdef __APPLE__
// cl.hpp is not included on OS X, have to rely on a custom cl.hpp file
// provided by viennacl
//#include <OpenCL/cl.hpp>
#include "cl.hpp"
#else
#include <CL/cl.hpp>
#endif

#define OPENCL_CHECK_ERR( err, name ) do {                                                                          \
   cl_int rerr = err;                                                                                                      \
   if( rerr != CL_SUCCESS ) {                                                                                        \
       std::cerr << "OpenCL error (" << rerr << " " << std::hex << rerr << ") in file '" << __FILE__ << " in line " << __LINE__ << ": " << name << std::endl;    \
       exit(EXIT_FAILURE);                                                                                          \
   } } while (0)

class crbCLContext
{
    private:
        /* Is the context created ? */
        bool initialized_;
        /* base cl context */
        cl::Context * context_;

        /* list of devices used by this context */
        std::vector< cl::Device > deviceList_;
        /* queues for each device */
        std::vector<cl::CommandQueue *> queue_;
        /* available programs */
        std::map<std::string, cl::Program *> program_;
        /* keep data about already allocated memory segments */
        std::map<std::string, cl::Buffer *> mempool_;
        /* ViennaCL does not provide a way to know if it has already been initialized */
        bool linkedToVCL_;
        
    public:
        crbCLContext()
        {
            initialized_ = false;
            context_ = NULL;
            linkedToVCL_ = false;
        }

        /* Create the OpenCL context for a device type and required device count */
        int init(cl_device_type type, int count)
        {
            cl_int err;
            /* list of available platforms on current server */
            std::vector< cl::Platform > platformList;
            std::vector< cl::Device > allList;
            std::vector< cl::Device > tList;

            if(!initialized_)
            {
                /* TODO cleaning stuff if we realloc context */

                /* Get available OpenCL platforms */
                cl::Platform::get(&platformList);
                LOG( INFO ) << "[CRB::fixedPointPrimalCL] Platform count is: " << platformList.size() << std::endl;
                err = platformList.size()!=0 ? CL_SUCCESS : -1;
                OPENCL_CHECK_ERR(err, "cl::Platform::get: No OpenCL Platform available");

                /* Gather available devices on the current node */
                for(size_t k = 0; k < platformList.size(); k++)
                {
                    std::string platformVendor, platformName;
                    platformList[k].getInfo((cl_platform_info)CL_PLATFORM_VENDOR, &platformVendor);
                    platformList[k].getInfo((cl_platform_info)CL_PLATFORM_NAME, &platformName);
                    //std::cout << "Platform " << k << " (" << platformName << ") is by: " << platformVendor << "\n";

                    try {
                        platformList[k].getDevices(CL_DEVICE_TYPE_ALL, &tList);
                        allList.insert(allList.end(), tList.begin(), tList.end());
                        tList.clear();
                    }
                    catch(cl::Error error) {
                        std::cout << error.what() << "(" << error.err() << ")" << std::endl;
                    }
                }

                /* Check for device availability */
                cl_device_type devType = false;
                cl_bool devAvail = false;
                std::string dname;
                for(int devID = 0; devID < allList.size(); devID++)
                {
                    try {
                        allList[devID].getInfo(CL_DEVICE_TYPE, &devType);
                        allList[devID].getInfo(CL_DEVICE_AVAILABLE, &devAvail);
                    }
                    catch(cl::Error error) {
                        std::cout << error.what() << "(" << error.err() << ")" << std::endl;
                    }

                    /* if we found a device available and of the required type */
                    /* we create a queue for it */
                    if(devType == type && devAvail)
                    {
                        deviceList_.push_back(allList[devID]);
                    }

                    /* if the number of queues is equal to the required number of device */
                    /* then we break the loop */
                    if(deviceList_.size() == count)
                    {
                        break;
                    }
                }

                /* create the cl context */
                context_ = new cl::Context(deviceList_, NULL, NULL, NULL, &err);

                /* create the corresponding queues */
                for(int i = 0; i < deviceList_.size(); i++)
                {
                    queue_.push_back(new cl::CommandQueue(*context_, deviceList_[i], 0, &err));

                    deviceList_[i].getInfo(CL_DEVICE_TYPE, &devType);
                    deviceList_[i].getInfo(CL_DEVICE_NAME, &dname);
                    std::cout << "Using device " << i << ":" << std::endl;
                    std::cout << "- Name: " << dname << std::endl;
                    std::cout << "- Type: " << devType << std::endl;

                    /*
                    cl_uint maxCU;
                    deviceList_[i].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &maxCU);
                    std::cout << "- CL_DEVICE_MAX_COMPUTE_UNITS: " << maxCU << std::endl;

                    cl_uint maxWIDim;
                    deviceList_[i].getInfo(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &maxWIDim);
                    std::cout << "- CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: " << maxWIDim << std::endl;

                    size_t * mwis = new size_t[maxWIDim];
                    deviceList_[i].getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, mwis);
                    for(int j = 0; j < maxWIDim; j++)
                    {
                        std::cout << " * Max work items[" << j << "]: " << mwis[j] << std::endl;
                    }
                    delete[] mwis;

                    size_t maxWGsize;
                    deviceList_[i].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &maxWGsize);
                    std::cout << "- CL_DEVICE_MAX_WORK_GROUP_SIZE: " << maxWGsize << std::endl;
                    */
                }

                initialized_ = true;
            }

            return queue_.size();
        }

        //void setContext(int devID, cl::Context * ctx){ if(!context_){ devID_ = devID; context_ = ctx;} }
        //int getDeviceID(){ return devID_; }
        cl::Context * getContext(){ return context_; }

        cl::Device * getDevice(int index){ if(index < deviceList_.size()) { return &(deviceList_[index]); } return NULL; }

        //void setCommandQueue(cl::CommandQueue * queue){ if(!queue_){ queue_ = queue;} }
        cl::CommandQueue * getCommandQueue(int index){ if(index < queue_.size()) { return queue_[index]; } return NULL; }

        cl::Program * getProgram(std::string name)
        {
            if(program_.find(name) != program_.end())
            {
                return program_[name];
            }

            return NULL;
        }

        bool addProgram(std::string name, const char * src, size_t size)
        {
            cl_int err;

            /* The program already exists ? */
            if(program_.find("name") != program_.end())
            {
                return false;
            }

            /* declare sources */
            cl::Program::Sources source(1, std::make_pair(src, size));

            cl::Program * program = new cl::Program(*context_, source, &err);
            OPENCL_CHECK_ERR(err, "Could not init program");
            err = program->build(deviceList_);
            OPENCL_CHECK_ERR(err, "Could not build kernel");
            if(err != CL_SUCCESS)
            {
                cl_build_status status;
                program->getBuildInfo(deviceList_.front(), CL_PROGRAM_BUILD_STATUS, &status);

                std::string log;
                program->getBuildInfo(deviceList_.front(), CL_PROGRAM_BUILD_LOG, &log);
                std::cout << log  << std::endl;

                delete program;
                return false;
            }
            else
            {
                program_[name] = program;
            }

            return true;
        }

        cl::Buffer * getBuffer(std::string name, cl_mem_flags flags, size_t size, void* host_ptr = NULL, cl_int* err = NULL)
        {
            cl::Buffer * B = NULL;

            /* If the error container is specified */
            /* we set it by default to CL_SUCCESS */
            /* for the cases that don't use CL functions to work */
            if(err) {
                *err = CL_SUCCESS;
            }

            std::map<std::string, cl::Buffer *>::iterator it;
            /* Try to recover a previous named  buffer form the memory pool */
            it = mempool_.find(name);
            if(it != mempool_.end()) {
                //std::cout << name << ": Found existing buffer" << std::endl;
                /* if the buffer has the same characteristics */
                /* we can directly use it */
                cl_mem_flags pflags;
                it->second->getInfo(CL_MEM_FLAGS, &pflags);
                size_t psize;
                it->second->getInfo(CL_MEM_SIZE, &psize);

                /* same buffer */
                if(size == psize && flags == pflags) {
                    //std::cout << name << ": Using the same buffer" << std::endl;
                    B = it->second;
                }
                /* different buffer => delete and recreate */
                else {
                    //std::cout << name << ": Reallocating new buffer (" << size << "/" << psize << ", " << flags << "/" << pflags << ")" << std::endl;
                    delete it->second;
                    mempool_[name] = B = new cl::Buffer(*context_, flags, size, host_ptr, err);
                }
            }
            /* otherwise we alloc a new buffer */
            else
            {
                //std::cout << name << ": Allocating new buffer" << std::endl;
                mempool_[name] = B = new cl::Buffer(*context_, flags, size, host_ptr, err);
            }

            return B;
        }

        bool isLinkedToVCL(){ return linkedToVCL_; }
        void setLinkedToVCL(){ linkedToVCL_ = true; }

        void clean()
        {
            /* free buffers */
            std::map<std::string, cl::Buffer *>::iterator itm;
            for(itm = mempool_.begin(); itm != mempool_.end(); itm++)
            {
                if(itm->second)
                {
                    delete itm->second;
                    itm->second = NULL;
                }
            }
            mempool_.clear();

            /* clean programs created */
            std::map<std::string, cl::Program *>::iterator itp;
            for(itp = program_.begin(); itp != program_.end(); itp++)
            {
                if(itp->second)
                {
                    delete itp->second;
                    itp->second = NULL;
                }
            }
            program_.clear();

            /* clean device queues created */
            for(int i = 0; i < queue_.size(); i++)
            {
                delete queue_[i];
            }
            queue_.clear();

            delete context_;
            context_ = NULL;
        }

        ~crbCLContext()
        {
            this->clean();
        }
};


#endif
