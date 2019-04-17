// ###############################################################################################
// planeResolSq: plane normal tolerance
// misorResol: misorientation tolerance, completely depends on misorientation angle. 


// ###############################################################################################
// Get symmetricall equivalent grain_normals
for(int j = 0; j < nsym; j++)
{
  // rotate g1 by symOp
  m_OrientationOps[cryst]->getMatSymOp(j, sym1);
  MatrixMath::Multiply3x3with3x3(sym1, g1, g1s);
  // get the crystal directions along the triangle normals
  MatrixMath::Multiply3x3with3x1(g1s, normal_lab, normal_grain1);

  for(int k = 0; k < nsym; k++)
  {
    // calculate the symmetric misorienation
    m_OrientationOps[cryst]->getMatSymOp(k, sym2);
    // rotate g2 by symOp
    MatrixMath::Multiply3x3with3x3(sym2, g2, g2s);
    // transpose rotated g2
    MatrixMath::Transpose3x3(g2s, g2sT);
    // calculate delta g
    MatrixMath::Multiply3x3with3x3(g1s, g2sT, dg); // dg -- the misorientation between adjacent grains
    MatrixMath::Transpose3x3(dg, dgT);

    for(int transpose = 0; transpose <= 1; transpose++)
    {
      // check if dg is close to gFix
      if(transpose == 0)
      {
        MatrixMath::Multiply3x3with3x3(dg, gFixedT, diffFromFixed);
      }
      else
      {
        MatrixMath::Multiply3x3with3x3(dgT, gFixedT, diffFromFixed);
      }

      float diffAngle = acosf((diffFromFixed[0][0] + diffFromFixed[1][1] + diffFromFixed[2][2] - 1.0f) * 0.5f);

      if(diffAngle < m_misorResol)
      {
        MatrixMath::Multiply3x3with3x1(dgT, normal_grain1, normal_grain2); // minus sign before normal_grain2 will be added later

        if(transpose == 0)
        {
          (*selectedTris).push_back(TriAreaAndNormals(m_FaceAreas[triIdx], normal_grain1[0], normal_grain1[1], normal_grain1[2], -normal_grain2[0], -normal_grain2[1], -normal_grain2[2]));
        }
        else
        {
          (*selectedTris).push_back(TriAreaAndNormals(m_FaceAreas[triIdx], -normal_grain2[0], -normal_grain2[1], -normal_grain2[2], normal_grain1[0], normal_grain1[1], normal_grain1[2]));
        }


// ###############################################################################################
// Compare grain_normals to fixed_normals
  virtual ~ProbeDistrib() = default;

  void probe(size_t start, size_t end) const
  {
    for(size_t ptIdx = start; ptIdx < end; ptIdx++)
    {
      float fixedNormal1[3] = {samplPtsX.at(ptIdx), samplPtsY.at(ptIdx), samplPtsZ.at(ptIdx)};
      float fixedNormal2[3] = {0.0f, 0.0f, 0.0f};
      MatrixMath::Multiply3x3with3x1(gFixedT, fixedNormal1, fixedNormal2);

      for(int triRepresIdx = 0; triRepresIdx < static_cast<int>(selectedTris.size()); triRepresIdx++)
      {
        for(int inversion = 0; inversion <= 1; inversion++)
        {
          float sign = 1.0f;
          if(inversion == 1)
          {
            sign = -1.0f;
          }

          float theta1 = acosf(sign * (selectedTris[triRepresIdx].normal_grain1_x * fixedNormal1[0] + selectedTris[triRepresIdx].normal_grain1_y * fixedNormal1[1] +
                                       selectedTris[triRepresIdx].normal_grain1_z * fixedNormal1[2]));

          float theta2 = acosf(-sign * (selectedTris[triRepresIdx].normal_grain2_x * fixedNormal2[0] + selectedTris[triRepresIdx].normal_grain2_y * fixedNormal2[1] +
                                        selectedTris[triRepresIdx].normal_grain2_z * fixedNormal2[2]));

          float distSq = 0.5f * (theta1 * theta1 + theta2 * theta2);

          if(distSq < planeResolSq)
          {
            (*distribValues)[ptIdx] += selectedTris[triRepresIdx].area;
          }
        }
      }
      (*errorValues)[ptIdx] = sqrt((*distribValues)[ptIdx] / totalFaceArea / double(numDistinctGBs)) / ballVolume;

      (*distribValues)[ptIdx] /= totalFaceArea;
      (*distribValues)[ptIdx] /= ballVolume;

// ----- ballVolume -----
const float FindGBCDMetricBased::k_ResolutionChoices[FindGBCDMetricBased::k_NumberResolutionChoices][2] = {{3.0f, 7.0f}, {5.0f, 5.0f}, {5.0f, 7.0f}, {5.0f, 8.0f},
                                                                                                           {6.0f, 7.0f}, {7.0f, 7.0f}, {8.0f, 8.0f}}; // { for misorient., for planes }
const double FindGBCDMetricBased::k_BallVolumesM3M[FindGBCDMetricBased::k_NumberResolutionChoices] = {0.0000641361, 0.000139158, 0.000287439, 0.00038019, 0.000484151, 0.000747069, 0.00145491};


// ###############################################################################################
// Plane normal to Aenith angles
float zenith = acosf(samplPtsZ.at(ptIdx));
float azimuth = atan2f(samplPtsY.at(ptIdx), samplPtsX.at(ptIdx));
fprintf(fDist, "%.2f %.2f %.4f\n", azimuthDeg, 90.0f - zenithDeg, distribValues[ptIdx]);



// ###############################################################################################
// Get misorientation matrix 
float gFixed[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
float gFixedT[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};

{
  float gFixedAngle = static_cast<float>(m_MisorientationRotation.angle * SIMPLib::Constants::k_PiOver180);
  float gFixedAxis[3] = {m_MisorientationRotation.h, m_MisorientationRotation.k, m_MisorientationRotation.l};
  MatrixMath::Normalize3x1(gFixedAxis);
  FOrientArrayType om(9);
  FOrientTransformsType::ax2om(FOrientArrayType(gFixedAxis[0], gFixedAxis[1], gFixedAxis[2], gFixedAngle), om);
  om.toGMatrix(gFixed);
}

MatrixMath::Transpose3x3(gFixed, gFixedT);


 // ###############################################################################################
// generate uniformly sampled points, see http://www.softimageblog.com/archives/115
int numSamplPts_WholeSph = 2 * m_NumSamplPts; // here we generate points on the whole sphere
QVector<float> samplPtsX(0);
QVector<float> samplPtsY(0);
QVector<float> samplPtsZ(0);

float _inc = 2.3999632f; // = pi * (3 - sqrt(5))
float _off = 2.0f / float(numSamplPts_WholeSph);

for(int ptIdx_WholeSph = 0; ptIdx_WholeSph < numSamplPts_WholeSph; ptIdx_WholeSph++)
{
  if(getCancel())
  {
    return;
  }

  float _y = (float(ptIdx_WholeSph) * _off) - 1.0f + (0.5f * _off);
  float _r = sqrtf(fmaxf(1.0f - _y * _y, 0.0f));
  float _phi = float(ptIdx_WholeSph) * _inc;

  float z = sinf(_phi) * _r;

  if(z > 0.0f)
  {
    samplPtsX.push_back(cosf(_phi) * _r);
    samplPtsY.push_back(_y);
    samplPtsZ.push_back(z);
  }
}
// Add points at the equator for better performance of some plotting tools
for(double phi = 0.0; phi <= SIMPLib::Constants::k_2Pi; phi += m_planeResol)
{
  samplPtsX.push_back(cosf(static_cast<float>(phi)));
  samplPtsY.push_back(sinf(static_cast<float>(phi)));
  samplPtsZ.push_back(0.0f);
}

