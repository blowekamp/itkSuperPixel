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
#ifndef itkSLICImageFilter_h
#define itkSLICImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkBarrier.h"


namespace itk
{

template< typename TInputImage, typename TOutputImage, typename TDistancePixel = float>
class SLICImageFilter:
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(SLICImageFilter);

  /** Standard class type aliases. */
  using Self = SLICImageFilter;
  using Superclass = ImageToImageFilter< TInputImage, TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SLICImageFilter, ImageToImageFilter);

  /** ImageDimension constants */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;


  /** Image type information. */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputImageType = TOutputImage;
  using DistanceType = TDistancePixel;
  using DistanceImageType = Image<DistanceType, ImageDimension>;

  using IndexType = typename InputImageType::IndexType;
  using PointType = typename InputImageType::PointType;

  // assume variable length vector right now
  using ClusterComponentType = double;
  using ClusterType = vnl_vector_ref<ClusterComponentType>;

  using OutputImageRegionType = typename OutputImageType::RegionType;

  using SuperGridSizeType = FixedArray< unsigned int, ImageDimension >;

  itkSetMacro( SpatialProximityWeight, double );
  itkGetConstMacro( SpatialProximityWeight, double );

  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetConstMacro( MaximumNumberOfIterations, unsigned int );

  itkSetMacro(SuperGridSize, SuperGridSizeType);
  void SetSuperGridSize(unsigned int factor);
  void SetSuperGridSize(unsigned int i, unsigned int factor);

  itkSetMacro(EnforceConnectivity, bool);
  itkGetMacro(EnforceConnectivity, bool);
  itkBooleanMacro(EnforceConnectivity);


  /* Get the current average cluster residual.
   *
   * After each iteration the residual is computed as the distance
   * between the current clusters and the previous. This is averaged
   * so that the value is independent of the number of clusters.
   */
  itkGetConstMacro( AverageResidual, double );

protected:
  SLICImageFilter();
  ~SLICImageFilter();

  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Generate full output and require full input */
  void EnlargeOutputRequestedRegion(DataObject *output) override;

  void BeforeThreadedGenerateData() override;

  void ThreadedUpdateDistanceAndLabel(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

  void ThreadedUpdateClusters(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

  void ThreadedPerturbClusters(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

  void ThreadedConnectivity(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;



  void AfterThreadedGenerateData() override;

  DistanceType Distance(const ClusterType &cluster1, const ClusterType &cluster2);

  DistanceType Distance(const ClusterType &cluster, const InputPixelType &v, const PointType &pt);

private:

  SuperGridSizeType m_SuperGridSize;
  unsigned int      m_MaximumNumberOfIterations;
  double            m_SpatialProximityWeight;

  FixedArray<double,ImageDimension> m_DistanceScales;
  std::vector<ClusterComponentType> m_Clusters;
  std::vector<ClusterComponentType> m_OldClusters;

  struct UpdateCluster
  {
    size_t count;
    vnl_vector<ClusterComponentType> cluster;
  };

  typedef std::map<size_t, UpdateCluster> UpdateClusterMap;

  std::vector<UpdateClusterMap> m_UpdateClusterPerThread;

  typename Barrier::Pointer           m_Barrier;
  typename DistanceImageType::Pointer m_DistanceImage;

  bool m_EnforceConnectivity;

  double m_AverageResidual;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSLICImageFilter.hxx"
#endif

#endif //itkSLICImageFilter_h
