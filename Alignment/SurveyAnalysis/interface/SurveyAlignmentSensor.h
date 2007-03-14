#ifndef Alignment_SurveyAnalysis_SurveyAlignmentSensor_h
#define Alignment_SurveyAnalysis_SurveyAlignmentSensor_h

/** \class SurveyAlignmentSensor
 *
 *  Survey alignment using sensor residual.
 *
 *  The residual (dRx, dRy, dRz, dWx, dWy, dWz) is found for each sensor.
 *  The sensor shifted by this amount during each iteration.
 *
 *  $Date: 2007/02/14 $
 *  $Revision: 1 $
 *  \author Chung Khim Lae
 */

#include "Alignment/SurveyAnalysis/interface/SurveyAlignment.h"

class SurveyAlignmentSensor:
  public SurveyAlignment
{
  public:

  /// Constructor to set the sensors in base class.
  SurveyAlignmentSensor(
			const std::vector<Alignable*>& sensors
			);

  protected:

  /// Find the alignment parameters for all sensors.
  virtual void findAlignPars();
};

#endif
