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
#ifndef itkSLICImageFilter_hxx
#define itkSLICImageFilter_hxx

#include "itkSLICImageFilter.h"

namespace itk
{

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::SLICImageFilter()
  : m_MaximumNumberOfIterations( (ImageDimension > 2) ? 5 : 10),
    m_SpatialProximityWeight( 10.0 ),
    m_Barrier(Barrier::New())
{
  m_SuperGridSize.Fill(50);
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::~SLICImageFilter()
{
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::SetSuperGridSize(unsigned int factor)
{
  unsigned int i;
  for (i = 0; i < ImageDimension; ++i)
    {
    if (factor != m_SuperGridSize[i])
      {
      break;
      }
    }
  if ( i < ImageDimension )
    {
    this->Modified();
    m_SuperGridSize.Fill(factor);
    }
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::SetSuperGridSize(unsigned int i, unsigned int factor)
{
  if (m_SuperGridSize[i] == factor)
    {
    return;
    }

  this->Modified();
  m_SuperGridSize[i] = factor;
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "SuperGridSize: " << m_SuperGridSize << std::endl;
  os << indent << "MaximumNumberOfIterations: " << m_MaximumNumberOfIterations << std::endl;
  os << indent << "SpatialProximityWeight: " << m_SpatialProximityWeight << std::endl;
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::BeforeThreadedGenerateData() ITK_OVERRIDE
{
  itkDebugMacro("Starting BeforeThreadedGenerateData");


  ThreadIdType numberOfThreads = this->GetNumberOfThreads();

  if ( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 )
    {
    numberOfThreads = vnl_math_min(
      this->GetNumberOfThreads(), itk::MultiThreader::GetGlobalMaximumNumberOfThreads() );
    }

  // number of threads can be constrained by the region size, so call the
  // SplitRequestedRegion to get the real number of threads which will be used
  typename TOutputImage::RegionType splitRegion;  // dummy region - just to call
  // the following method

  numberOfThreads = this->SplitRequestedRegion(0, numberOfThreads, splitRegion);

  m_Barrier->Initialize(numberOfThreads);

  const InputImageType *inputImage = this->GetInput();

  itkDebugMacro("Shinking Starting")
    typename InputImageType::Pointer shrunkImage;
  {
  // todo disconnect input from pipeline
  typedef itk::ShrinkImageFilter<InputImageType, InputImageType> ShrinkImageFilterType;
  typename ShrinkImageFilterType::Pointer shrinker = ShrinkImageFilterType::New();
  shrinker->SetInput(inputImage);
  shrinker->SetShrinkFactors(m_SuperGridSize);
  shrinker->UpdateLargestPossibleRegion();

  shrunkImage = shrinker->GetOutput();
  }
  itkDebugMacro("Shinking Completed")

  const typename InputImageType::RegionType region = inputImage->GetBufferedRegion();
  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;
  const size_t numberOfClusters = shrunkImage->GetBufferedRegion().GetNumberOfPixels();


  // allocate array of scalars
  m_Clusters.resize(numberOfClusters*numberOfClusterComponents);

  typedef ImageScanlineConstIterator< InputImageType > InputConstIteratorType;


  InputConstIteratorType it(shrunkImage, shrunkImage->GetLargestPossibleRegion());

  // Initialize cluster centers
  size_t cnt = 0;
  while(!it.IsAtEnd())
    {
    const size_t         ln =  shrunkImage->GetLargestPossibleRegion().GetSize(0);
    for (unsigned x = 0; x < ln; ++x)
      {
      // construct vector as reference to the scalar array
      ClusterType cluster( numberOfClusterComponents, &m_Clusters[cnt*numberOfClusterComponents] );

      for(unsigned int i = 0; i < numberOfComponents; ++i)
        {
        cluster[i] = it.Get()[i];
        }
      const IndexType & idx = it.GetIndex();
      typename InputImageType::PointType pt;
      shrunkImage->TransformIndexToPhysicalPoint(idx, pt);
      for(unsigned int i = 0; i < ImageDimension; ++i)
        {
        cluster[numberOfComponents+i] = pt[i];
        }
      ++it;
      ++cnt;
      }
    it.NextLine();
    }
  itkDebugMacro("Initial Clustering Completed");

  shrunkImage = ITK_NULLPTR;

  // TODO: Move cluster center to lowest gradient position in a 3x
  // neighborhood


  m_DistanceImage = DistanceImageType::New();
  m_DistanceImage->CopyInformation(inputImage);
  m_DistanceImage->SetBufferedRegion( region );
  m_DistanceImage->Allocate();

  for (unsigned int i = 0; i < ImageDimension; ++i)
    {
    m_DistanceScales[i] = m_SpatialProximityWeight/m_SuperGridSize[i];
    }


  // deep copy to ensure memory is allocated
  std::vector<ClusterComponentType>(m_Clusters.begin(), m_Clusters.end()).swap(m_OldClusters);

  this->Superclass::BeforeThreadedGenerateData();
}



template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::ThreadedUpdateDistanceAndLabel(const OutputImageRegionType & outputRegionForThread, ThreadIdType  itkNotUsed(threadId))
{
  typedef ImageScanlineConstIterator< InputImageType > InputConstIteratorType;
  typedef ImageScanlineIterator< DistanceImageType >   DistanceIteratorType;

  const InputImageType *inputImage = this->GetInput();
  OutputImageType *outputImage = this->GetOutput();
  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;

  typename InputImageType::SizeType searchRadius;
  for (unsigned int i = 0; i < ImageDimension; ++i)
    {
    searchRadius[i] = m_SuperGridSize[i];
    }

  for (size_t i = 0; i*numberOfClusterComponents < m_Clusters.size(); ++i)
    {
    ClusterType cluster(numberOfClusterComponents, &m_Clusters[i*numberOfClusterComponents]);
    typename InputImageType::RegionType localRegion;
    typename InputImageType::PointType pt;
    IndexType idx;

    for (unsigned int d = 0; d < ImageDimension; ++d)
      {
      pt[d] = cluster[numberOfComponents+d];
      }
    //std::cout << "Cluster " << i << "@" << pt <<": " << cluster << std::endl;
    inputImage->TransformPhysicalPointToIndex(pt, idx);

    localRegion.SetIndex(idx);
    localRegion.GetModifiableSize().Fill(1u);
    localRegion.PadByRadius(searchRadius);
    if (!localRegion.Crop(outputRegionForThread))
      {
      continue;
      }


    const size_t         ln =  localRegion.GetSize(0);

    InputConstIteratorType inputIter(inputImage, localRegion);
    DistanceIteratorType   distanceIter(m_DistanceImage, localRegion);


    while ( !inputIter.IsAtEnd() )
      {
      for( size_t x = 0; x < ln; ++x )
        {
        const IndexType &currentIdx = inputIter.GetIndex();

        inputImage->TransformIndexToPhysicalPoint(currentIdx, pt);
        const double distance = this->Distance(cluster,
                                               inputIter.Get(),
                                               pt);
        if (distance < distanceIter.Get() )
          {
          distanceIter.Set(distance);
          outputImage->SetPixel(currentIdx, i);
          }

        ++distanceIter;
        ++inputIter;
        }
      inputIter.NextLine();
      distanceIter.NextLine();
      }

    // for neighborhood iterator size S
    }

}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId)
{
  typedef typename InputImageType::PixelType InputPixelType;

  const InputImageType *inputImage = this->GetInput();
  OutputImageType *outputImage = this->GetOutput();

  const typename InputImageType::RegionType region = inputImage->GetBufferedRegion();
  const unsigned int numberOfComponents = inputImage->GetNumberOfComponentsPerPixel();
  const unsigned int numberOfClusterComponents = numberOfComponents+ImageDimension;

  typedef ImageScanlineConstIterator< InputImageType > InputConstIteratorType;
  typedef ImageScanlineIterator< OutputImageType >     OutputIteratorType;


  itkDebugMacro("Entering Main Loop");
  for(unsigned int loopCnt = 0;  loopCnt<m_MaximumNumberOfIterations; ++loopCnt)
    {
    itkDebugMacro("Iteration :" << loopCnt);

    if (threadId == 0)
      {
      m_DistanceImage->FillBuffer(NumericTraits<typename DistanceImageType::PixelType>::max());
      }
    m_Barrier->Wait();

    ThreadedUpdateDistanceAndLabel(outputRegionForThread,threadId);

    m_Barrier->Wait();

    if (threadId==0)
      {
      // clear
      swap(m_Clusters, m_OldClusters);
      std::fill(m_Clusters.begin(), m_Clusters.end(), 0.0);


      std::vector<size_t> clusterCount(m_Clusters.size()/numberOfClusterComponents, 0);
      itkDebugMacro("Estimating Centers");
      // calculate new centers
      OutputIteratorType itOut = OutputIteratorType(outputImage, region);
      InputConstIteratorType itIn = InputConstIteratorType(inputImage, region);
      while(!itOut.IsAtEnd()&& threadId==0)
        {
        const size_t         ln =  region.GetSize(0);
        for (unsigned x = 0; x < ln; ++x)
          {
          const IndexType &idx = itOut.GetIndex();
          const InputPixelType &v = itIn.Get();
          const typename OutputImageType::PixelType l = itOut.Get();

          ClusterType cluster(numberOfClusterComponents, &m_Clusters[l*numberOfClusterComponents]);
          ++clusterCount[l];

          for(unsigned int i = 0; i < numberOfComponents; ++i)
            {
            cluster[i] += v[i];
            }

          typename InputImageType::PointType pt;
          inputImage->TransformIndexToPhysicalPoint(idx, pt);
          for(unsigned int i = 0; i < ImageDimension; ++i)
            {
            cluster[numberOfComponents+i] += pt[i];
            }

          ++itIn;
          ++itOut;
          }
        itIn.NextLine();
        itOut.NextLine();
        }

      // average, l1
      double l1Residual = 0.0;
      for (size_t i = 0; i*numberOfClusterComponents < m_Clusters.size(); ++i)
        {
        ClusterType cluster(numberOfClusterComponents,&m_Clusters[i*numberOfClusterComponents]);
        cluster /= clusterCount[i];

        ClusterType oldCluster(numberOfClusterComponents, &m_OldClusters[i*numberOfClusterComponents]);
        l1Residual += Distance(cluster,oldCluster);

        }

      std::cout << "L1 residual: " << std::sqrt(l1Residual) << std::endl;
      }
    // while error <= threshold
    }


}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
void
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::AfterThreadedGenerateData() ITK_OVERRIDE
{
  itkDebugMacro("Starting AfterThreadedGenerateData");


  // cleanup
  std::vector<ClusterComponentType>().swap(m_Clusters);
  std::vector<ClusterComponentType>().swap(m_OldClusters);
}


template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
typename SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>::DistanceType
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::Distance(const ClusterType &cluster1, const ClusterType &cluster2)
{
  const unsigned int s = cluster1.size();
  DistanceType d1 = 0.0;
  DistanceType d2 = 0.0;
  unsigned int i = 0;
  for (; i<s-ImageDimension; ++i)
    {
    const DistanceType d = (cluster1[i] - cluster2[i]);
    d1 += d*d;
    }
  //d1 = std::sqrt(d1);

  for (unsigned int j = 0; j < ImageDimension; ++j)
    {
    const DistanceType d = (cluster1[i] - cluster2[i]) * m_DistanceScales[j];
    d2 += d*d;
    ++i;
    }
  //d2 = std::sqrt(d2);
  return d1+d2;
}

template<typename TInputImage, typename TOutputImage, typename TDistancePixel>
typename SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>::DistanceType
SLICImageFilter<TInputImage, TOutputImage, TDistancePixel>
::Distance(const ClusterType &cluster, const InputPixelType &v, const PointType &pt)
{
  const unsigned int s = cluster.size();
  DistanceType d1 = 0.0;
  DistanceType d2 = 0.0;
  unsigned int i = 0;
  for (; i<s-ImageDimension; ++i)
    {
    const DistanceType d = (cluster[i] - v[i]);
    d1 += d*d;
    }
  //d1 = std::sqrt(d1);

  for (unsigned int j = 0; j < ImageDimension; ++j)
    {
    const DistanceType d = (cluster[i] - pt[j]) * m_DistanceScales[j];
    d2 += d*d;
    ++i;
    }
  //d2 = std::sqrt(d2);
  return d1+d2;
}



} // end namespace itk

#endif // itkSLICImageFilter_hxx
