/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkGTest.h"

#include "itkSLICImageFilter.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCommand.h"

namespace
{

class SLICFixture
  : public ::testing::Test
{
public:

  SLICFixture() {}
  ~SLICFixture() override {}

protected:

  template<unsigned int D, typename TPixelType = unsigned short>
  struct FixtureUtilities
  {
    static const unsigned int Dimension = D;

    using PixelType = TPixelType;
    using OutputPixelType = unsigned int;
    using InputImageType = itk::Image<PixelType, Dimension>;
    using OutputImageType = itk::Image<OutputPixelType, Dimension>;

    using FilterType = itk::SLICImageFilter<InputImageType, OutputImageType>;
  };

};
}


TEST_F(SLICFixture, SetGetPrint)
{
  using namespace itk::GTest::TypedefsAndConstructors::Dimension3;
  using Utils = FixtureUtilities<3>;

  auto filter = Utils::FilterType::New();
  filter->Print(std::cout);

  typename Utils::FilterType::ConstPointer constfilter = (const Utils::FilterType*)(filter.GetPointer());

  EXPECT_STREQ("SLICImageFilter", filter->GetNameOfClass());
  EXPECT_STREQ("ImageToImageFilter",  filter->Superclass::GetNameOfClass());

  Utils::FilterType:: SuperGridSizeType gridSize3(3);
  EXPECT_NO_THROW(filter->SetSuperGridSize(gridSize3));
  EXPECT_VECTOR_NEAR(gridSize3, filter->GetSuperGridSize(), 0);

  EXPECT_NO_THROW(filter->SetSuperGridSize(4));
  EXPECT_VECTOR_NEAR(Utils::FilterType:: SuperGridSizeType(4), filter->GetSuperGridSize(), 0);

  EXPECT_NO_THROW(filter->SetMaximumNumberOfIterations(6));
  EXPECT_EQ(6, filter->GetMaximumNumberOfIterations());

  EXPECT_NO_THROW(filter->SetSpatialProximityWeight(9.1));
  EXPECT_EQ(9.1, filter->GetSpatialProximityWeight());

  EXPECT_NO_THROW(filter->EnforceConnectivityOn());
  EXPECT_TRUE(filter->GetEnforceConnectivity());
  EXPECT_NO_THROW(filter->EnforceConnectivityOff());
  EXPECT_FALSE(filter->GetEnforceConnectivity());

  EXPECT_NO_THROW(filter->SetEnforceConnectivity(true));
  EXPECT_TRUE(filter->GetEnforceConnectivity());
}
