# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component( MY_CURENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file( READ "${MY_CURENT_DIR}/README" DOCUMENTATION )

# itk_module() defines the module dependencies in ITKExternalTemplate
# ITKExternalTemplate depends on ITKCommon
# The testing module in ITKExternalTemplate depends on ITKTestKernel
# and ITKMetaIO(besides ITKExternalTemplate and ITKCore)

# define the dependencies of the include module and the tests
itk_module(ITKExternalTemplate
  DEPENDS
    ITKCommon
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_ALL
)
