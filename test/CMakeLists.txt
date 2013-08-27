itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(${itk-module}Tests
  itkEmptyTest.cxx
)

CreateTestDriver(${itk-module}  "${${itk-module}-Test_LIBRARIES}" "${${itk-module}Tests}")

itk_add_test(NAME itkDeleteMeEmptyTest
  COMMAND ${itk-module}TestDriver itkEmptyTest "argument1" "..." )
