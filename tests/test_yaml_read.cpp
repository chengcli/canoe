#include <gtest/gtest.h>

// extern int add(int a, int b);

class YamlReadTests : public testing::Test {
 public:
  YamlReadTests() {}

  virtual ~YamlReadTests() {}

  virtual void SetUp() {}

  virtual void TearDown() {}

  int add(int a, int b) { return a + b; }
};

TEST_F(YamlReadTests, test_add) {
  int actual = add(5, 10);
  int expected = 15;
  EXPECT_EQ(expected, actual);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
