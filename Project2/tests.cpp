#include <unittest++/UnitTest++.h>

/*TEST(WillFail){
    CHECK(false);
}*/

TEST(Potential1){
    CHECK(potential_HO(5.0)==25.0);
}

TEST(Potential2){
    CHECK(potential_C(1.,1.)==2.);
}

int main()
{
    return UnitTest::RunAllTests();
}
