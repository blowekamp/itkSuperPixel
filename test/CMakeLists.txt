itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(${itk-module}Tests
  itkSLICImageFilterTest.cxx
)

CreateTestDriver(${itk-module}  "${${itk-module}-Test_LIBRARIES}" "${${itk-module}Tests}")

itk_add_test(NAME itkSLICImageFilterTest1
  COMMAND ${itk-module}TestDriver itkSLICImageFilterTest
  "${CMAKE_CURRENT_SOURCE_DIR}/Data/Input/VM1111Shrink-RGB.png"
  "super1.png" 25 )


itk_add_test(NAME itkSLICImageFilterTest2
  COMMAND ${itk-module}TestDriver itkSLICImageFilterTest
  "${CMAKE_CURRENT_SOURCE_DIR}/Data/Input/VM1111Shrink-LAB.mha"
  "super2.png" 25 )
