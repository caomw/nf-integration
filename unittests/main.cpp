#include "unittests.h"

using namespace std;

int main() {

    CLinearAlgebraTest test0;
    QTest::qExec(&test0);
    CDenseArrayTest test1;
    QTest::qExec(&test1);
    return 0;

}

