/*
 * slimtest.h
 *
 * Created on: Mar 6, 2018
 *      Author: Tian Bo
 *
 * Constructing test case
 *
 * This is to provide a framework for small project - it does not make sense
 * to add a heavy test framework (e.g. gtest) in such case.
 *
 * In this framework, only one header file (this one) is required. There are
 * only three macros for constructing a test case: TEST, TEST_ASSERT
 * and TEST_FAIL. The syntax is exemplified with the following code
 *
 * TEST
 * {
 *     // test code here
 *     int a = 3;
 *     int b = 4;
 *
 *     TEST_ASSERT((a + b) == 7); // check test condition
 *
 *     if ((a + b) != (b + a)) {
 *         TEST_FAIL(); // mark a failure
 *     }
 * }
 *
 *
 * The main function can be simply the following:
 *
 * int main()
 * {
 *     tests.testAll();
 *     return 0;
 * }
 *
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>

#define CONCATE1(X,Y) X##Y
#define CONCATE(X,Y)  CONCATE1(X,Y)


struct TestBase
{
    virtual void testBody() = 0;
    virtual size_t index() = 0;
    virtual int line() = 0;
};

class Tests
{

public:
    struct Exception
    {
        int lineNumber;
    };

    ~Tests()
    {
        for(auto p : m_tests) {
            delete p;
        }
    }

    void add(TestBase * t)
    {
        m_tests.push_back(t);
    }

    size_t size()
    {
        return m_tests.size();
    }

    void testAll()
    {
        std::cout << "test line result" << std::endl;
        for(auto p : m_tests) {
            std::cout << std::setw(4) << p->index() << " "
                    << std::setw(4) << p->line();

            try {
                p->testBody();
                std::cout << " ok";
            }
            catch (Exception & e) {
                std::cout << " fail at " << e.lineNumber;
            }
            catch (std::exception & e) {
                std::cout << " exception: " << e.what();
            }
            catch (...) {
                std::cout << " unknown exception";
            }

            std::cout << std::endl;
        }
    }

private:
    std::vector<TestBase *> m_tests;
};

Tests tests;

#define TEST struct CONCATE(Test, __LINE__) : public TestBase \
{ \
    void testBody() final; \
    size_t index() final { return _index;} \
    int line() final {return _line;} \
    static size_t _index; \
    static size_t add() \
    { \
        tests.add(new CONCATE(Test, __LINE__)); \
        return tests.size(); \
    } \
    int _line = __LINE__; \
}; \
size_t CONCATE(Test, __LINE__)::_index = CONCATE(Test, __LINE__)::add(); \
void CONCATE(Test, __LINE__)::testBody()

#define TEST_ASSERT(cond) if (cond) {} else \
        throw Tests::Exception({__LINE__})

#define TEST_FAIL() TEST_ASSERT(false)
