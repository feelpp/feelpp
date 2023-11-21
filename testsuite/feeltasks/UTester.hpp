///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under LGPL Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef UTESTER_HPP
#define UTESTER_HPP

#include <iostream>
#include <list>
#include <string>
#include <cstdio>
#include <atomic>

template < class TestClass >
class UTester {
    // Test function pointer
    typedef void (TestClass::*TestFunc)(void);

    /** Test descriptor */
    struct TestFuncDescriptor {
        TestFunc func;    //< Test adress
        std::string name; //< Test name
    };


    std::list< TestFuncDescriptor > tests; //< all tests

    std::atomic<int> totalTests; //< number of tests

    std::atomic<int> currentTest; //< current processing test in the run
    std::atomic<int> currentStep; //< current processing step in the run

    std::atomic<int> failedSteps; //< number of failed step in the current test
    std::atomic<int> failedTests; //< number of failed tests

protected:
    /** Constructor */
    UTester() {
        totalTests = 0;
    }

    virtual ~UTester(){}

    /** Callback before processing test */
    virtual void Before() {
    }

    /** Callback after processing test */
    virtual void After() {
    }

    /** Callback before each unit test */
    virtual void PreTest() {
    }

    /** Callback after each unit test */
    virtual void PostTest() {
    }

    /**
    * This function has to add tests
        * <code> AddTest(&MyTest::TestOne); </code>
    */
    virtual void SetTests() = 0;

    /**
    * Add a test without giving a name
    * @param inFunc test function address
    */
    void AddTest(TestFunc inFunc) {
        char buff[256];
        sprintf(buff, "Unnamed Test number %d", totalTests + 1);
        AddTest(inFunc, buff);
    }

    /**
    * Add a test with a name
    * @param inFunc test function address
    * @param inFuncName function name
    */
    void AddTest(TestFunc inFunc, const std::string& inFuncName) {
        ++totalTests;
        TestFuncDescriptor desc;
        desc.func = inFunc;
        desc.name = inFuncName;
        tests.push_back(desc);
    }

    void uassertIsTrue(const bool shouldBeTrue,
                       const std::string str, const int line, const std::string file) {
        ++currentStep;
        if (!shouldBeTrue) {
            std::cout << ">> Step " << currentStep << " Equality Failed: " << str << "\n";
            std::cout << ">>     In file :" << file << "\n";
            std::cout << ">>     At line :" << line << "\n";
            ++failedSteps;
        }
    }

    template < class ValueType >
    void uassertAreEqual(const ValueType v1, const ValueType v2,
                         const std::string str, const int line, const std::string file) {
        ++currentStep;
        if (!(v1 == v2)) {
            std::cout << ">> Step " << currentStep << " Equality Failed: " << str << "\n";
            std::cout << ">>     In file :" << file << "\n";
            std::cout << ">>     At line :" << line << "\n";
            std::cout << ">>     v1 :" << v1 << "\n";
            std::cout << ">>     v2 :" << v2 << "\n";
            ++failedSteps;
        }
    }

    template < class ValueType >
    void uassertAreDiff(const ValueType v1, const ValueType v2,
                        const std::string str, const int line, const std::string file) {
        ++currentStep;
        if (!(v1 != v2)) {
            std::cout << ">> Step " << currentStep << " Difference Failed: " << str << "\n";
            std::cout << ">>     In file :" << file << "\n";
            std::cout << ">>     At line :" << line << "\n";
            std::cout << ">>     v1 :" << v1 << "\n";
            std::cout << ">>     v2 :" << v2 << "\n";
            ++failedSteps;
        }
    }

public:
    /**
    * Processing the test
        * return application exit code (= nb of errors)
    */
    int Run() {
        tests.clear();
        // register tests
        SetTests();

        TestClass* const toTest = static_cast< TestClass* >(this);
        currentTest             = 0;
        failedTests             = 0;

        Before();

        // for each tests
        const typename std::list< TestFuncDescriptor >::const_iterator end = tests.end();
        for (typename std::list< TestFuncDescriptor >::iterator iter = tests.begin(); iter != end; ++iter) {
            currentStep = 0;
            failedSteps = 0;

            std::cout << "[Start] " << (*iter).name << "\n";

            PreTest();
            TestFunc ff = (*iter).func;
            (toTest->*ff)();
            PostTest();

            if (failedSteps) {
                std::cout << "[Finished] FAILED (" << failedSteps << "/" << currentStep << " steps failed)\n";
                ++failedTests;
            } else {
                std::cout << "[Finished] PASSED (" << currentStep << " steps)\n";
            }

            ++currentTest;
        }


        After();

        std::cout << "Test is over, " << (totalTests - failedTests) << " Passed, " << failedTests << " Failed\n";

        return failedTests;
    }
};

#define UASSERTETRUE(condition) \
    this->uassertIsTrue(condition, #condition, __LINE__, __FILE__);

#define UASSERTEEQUAL(v1, v2) \
    this->uassertAreEqual(v1, v2, #v1 " == " #v2, __LINE__, __FILE__);

#define UASSERTEDIFF(v1, v2) \
    this->uassertAreDiff(v1, v2, #v1 " != " #v2, __LINE__, __FILE__);

#define UASSERTETRUE_OFFSET(condition, offset) \
    this->uassertIsTrue(condition, #condition, __LINE__ + offset, __FILE__);

#define UASSERTEEQUAL_OFFSET(v1, v2, offset) \
    this->uassertAreEqual(v1, v2, #v1 " == " #v2, __LINE__ + offset, __FILE__);

#define UASSERTEDIFF_OFFSET(v1, v2, offset) \
    this->uassertAreDiff(v1, v2, #v1 " != " #v2, __LINE__ + offset, __FILE__);

#endif

#define TestClass(X)\
    int main(void){\
    X Controller;\
    return Controller.Run();\
    }
