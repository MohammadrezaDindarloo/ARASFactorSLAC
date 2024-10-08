// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/rot3.h>

namespace sym {

/**
 * This function was autogenerated from a symbolic function. Do not modify by hand.
 *
 * Symbolic function: IK_residual_func_cost3_wrt_pd_Nl6
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     encoder: Matrix41
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix43
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPdNl6(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 570

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (187)
  const Scalar _tmp0 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = -_tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp7 = 2 * _tmp0;
  const Scalar _tmp8 = _tmp6 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp10 = _tmp2 * _tmp9;
  const Scalar _tmp11 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp8;
  const Scalar _tmp12 = 2 * _tmp2 * _tmp6;
  const Scalar _tmp13 = _tmp0 * _tmp9;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp5;
  const Scalar _tmp17 = _tmp16 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp18 = Scalar(1.0) / (_tmp17);
  const Scalar _tmp19 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp19 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp21 = -_tmp20;
  const Scalar _tmp22 = _tmp2 * _tmp7;
  const Scalar _tmp23 = _tmp6 * _tmp9;
  const Scalar _tmp24 =
      -Scalar(0.010999999999999999) * _tmp22 + Scalar(0.010999999999999999) * _tmp23;
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = _tmp21 + _tmp26;
  const Scalar _tmp28 = _tmp27 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp29 = _tmp18 * _tmp28;
  const Scalar _tmp30 =
      std::sqrt(Scalar(std::pow(_tmp17, Scalar(2)) + std::pow(_tmp28, Scalar(2))));
  const Scalar _tmp31 = Scalar(1.0) / (_tmp30);
  const Scalar _tmp32 = _tmp18 * _tmp30;
  const Scalar _tmp33 = _tmp32 * (_tmp16 * _tmp28 * _tmp31 - _tmp17 * _tmp27 * _tmp31);
  const Scalar _tmp34 =
      -Scalar(0.010999999999999999) * _tmp19 - Scalar(0.010999999999999999) * _tmp3;
  const Scalar _tmp35 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp36 = -_tmp35;
  const Scalar _tmp37 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp38 = _tmp34 + _tmp36 - _tmp37;
  const Scalar _tmp39 = _tmp29 * _tmp38;
  const Scalar _tmp40 = _tmp24 + _tmp25;
  const Scalar _tmp41 = _tmp20 + _tmp40;
  const Scalar _tmp42 = _tmp41 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp43 = _tmp11 + _tmp14;
  const Scalar _tmp44 = _tmp4 + _tmp43;
  const Scalar _tmp45 = _tmp44 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp46 = std::pow(Scalar(std::pow(_tmp42, Scalar(2)) + std::pow(_tmp45, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp47 = _tmp42 * _tmp46;
  const Scalar _tmp48 = _tmp45 * _tmp46;
  const Scalar _tmp49 = Scalar(1.0) / (_tmp29 * _tmp48 - _tmp47);
  const Scalar _tmp50 = _tmp34 + _tmp37;
  const Scalar _tmp51 = _tmp35 + _tmp50;
  const Scalar _tmp52 = _tmp38 * _tmp48;
  const Scalar _tmp53 = _tmp49 * (-_tmp29 * _tmp52 + _tmp47 * _tmp51);
  const Scalar _tmp54 = Scalar(1.0) * _tmp16;
  const Scalar _tmp55 = Scalar(1.0) * _tmp27;
  const Scalar _tmp56 = (-_tmp44 + _tmp54) / (_tmp41 - _tmp55);
  const Scalar _tmp57 = _tmp49 * (-_tmp48 * _tmp51 + _tmp52);
  const Scalar _tmp58 = _tmp29 * _tmp57 - _tmp38 - _tmp56 * (_tmp29 * _tmp53 + _tmp39);
  const Scalar _tmp59 = _tmp43 + _tmp5;
  const Scalar _tmp60 = _tmp59 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp61 = std::pow(_tmp60, Scalar(2));
  const Scalar _tmp62 = _tmp20 + _tmp26;
  const Scalar _tmp63 = _tmp62 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp64 = std::pow(_tmp63, Scalar(2));
  const Scalar _tmp65 = _tmp61 + _tmp64;
  const Scalar _tmp66 = std::pow(_tmp65, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp67 = _tmp38 * _tmp66;
  const Scalar _tmp68 = _tmp60 * _tmp67;
  const Scalar _tmp69 = _tmp36 + _tmp50;
  const Scalar _tmp70 = _tmp66 * _tmp69;
  const Scalar _tmp71 = _tmp63 * _tmp66;
  const Scalar _tmp72 = _tmp29 * _tmp66;
  const Scalar _tmp73 = _tmp60 * _tmp72 - _tmp71;
  const Scalar _tmp74 = -_tmp56 * (-_tmp29 * _tmp68 - _tmp53 * _tmp73 + _tmp63 * _tmp70) -
                        _tmp57 * _tmp73 - _tmp60 * _tmp70 + _tmp68;
  const Scalar _tmp75 = Scalar(1.0) / (_tmp74);
  const Scalar _tmp76 = _tmp62 * _tmp66;
  const Scalar _tmp77 = _tmp49 * (_tmp33 * _tmp48 + _tmp41 * _tmp48 - _tmp44 * _tmp47);
  const Scalar _tmp78 = _tmp33 * _tmp66;
  const Scalar _tmp79 = -_tmp59 * _tmp71 + _tmp60 * _tmp76 + _tmp60 * _tmp78 - _tmp73 * _tmp77;
  const Scalar _tmp80 = _tmp75 * _tmp79;
  const Scalar _tmp81 = _tmp29 * _tmp77 - _tmp33 - _tmp58 * _tmp80;
  const Scalar _tmp82 = Scalar(1.0) / (_tmp79);
  const Scalar _tmp83 = _tmp74 * _tmp82;
  const Scalar _tmp84 = _tmp58 + _tmp81 * _tmp83;
  const Scalar _tmp85 = _tmp73 * _tmp75;
  const Scalar _tmp86 = -_tmp29 - _tmp84 * _tmp85;
  const Scalar _tmp87 = _tmp48 * _tmp49;
  const Scalar _tmp88 = _tmp66 * _tmp75;
  const Scalar _tmp89 = _tmp84 * _tmp88;
  const Scalar _tmp90 = _tmp15 + _tmp4;
  const Scalar _tmp91 = _tmp90 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp92 = _tmp21 + _tmp40;
  const Scalar _tmp93 = _tmp92 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp94 = std::pow(Scalar(std::pow(_tmp91, Scalar(2)) + std::pow(_tmp93, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp95 = _tmp91 * _tmp94;
  const Scalar _tmp96 = _tmp32 * fh1;
  const Scalar _tmp97 = _tmp95 * _tmp96;
  const Scalar _tmp98 = Scalar(1.0) * _tmp53 * _tmp56 - Scalar(1.0) * _tmp57;
  const Scalar _tmp99 = -Scalar(1.0) * _tmp77 - _tmp80 * _tmp98;
  const Scalar _tmp100 = _tmp83 * _tmp99 + _tmp98;
  const Scalar _tmp101 = _tmp100 * _tmp88;
  const Scalar _tmp102 = -_tmp100 * _tmp85 + Scalar(1.0);
  const Scalar _tmp103 = _tmp93 * _tmp94;
  const Scalar _tmp104 = _tmp103 * _tmp96;
  const Scalar _tmp105 = _tmp54 + _tmp55 * _tmp56;
  const Scalar _tmp106 = 0;
  const Scalar _tmp107 = _tmp106 * _tmp88;
  const Scalar _tmp108 = _tmp106 * _tmp87;
  const Scalar _tmp109 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp110 = _tmp109 * _tmp32;
  const Scalar _tmp111 = Scalar(1.0) * _tmp82;
  const Scalar _tmp112 = _tmp111 * _tmp66;
  const Scalar _tmp113 = _tmp111 * _tmp87;
  const Scalar _tmp114 = fh1 * (_tmp103 * _tmp90 - _tmp92 * _tmp95);
  const Scalar _tmp115 = _tmp114 * _tmp32;
  const Scalar _tmp116 = std::exp(_tmp104 * (_tmp101 * _tmp60 + _tmp102 * _tmp87) +
                                  _tmp110 * (_tmp107 * _tmp60 - _tmp108 * _tmp85) +
                                  _tmp115 * (_tmp112 * _tmp60 - _tmp113 * _tmp73) +
                                  _tmp97 * (_tmp60 * _tmp89 + _tmp86 * _tmp87 + Scalar(1.0)));
  const Scalar _tmp117 = _tmp60 * _tmp66;
  const Scalar _tmp118 = std::pow(_tmp65, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp119 = _tmp118 * _tmp61;
  const Scalar _tmp120 = _tmp118 * _tmp60 * _tmp63;
  const Scalar _tmp121 = _tmp119 * _tmp29 - _tmp120 - _tmp72;
  const Scalar _tmp122 =
      _tmp119 * _tmp33 + _tmp119 * _tmp62 - _tmp120 * _tmp59 - _tmp121 * _tmp77 - _tmp76 - _tmp78;
  const Scalar _tmp123 = std::pow(_tmp79, Scalar(-2));
  const Scalar _tmp124 = Scalar(1.0) * _tmp123;
  const Scalar _tmp125 = _tmp122 * _tmp124;
  const Scalar _tmp126 = _tmp73 * _tmp87;
  const Scalar _tmp127 = _tmp119 * _tmp75;
  const Scalar _tmp128 = _tmp119 * _tmp38;
  const Scalar _tmp129 = _tmp120 * _tmp69;
  const Scalar _tmp130 =
      -_tmp119 * _tmp69 - _tmp121 * _tmp57 + _tmp128 -
      _tmp56 * (-_tmp121 * _tmp53 - _tmp128 * _tmp29 + _tmp129 + _tmp29 * _tmp67) - _tmp67 + _tmp70;
  const Scalar _tmp131 = std::pow(_tmp74, Scalar(-2));
  const Scalar _tmp132 = _tmp130 * _tmp131;
  const Scalar _tmp133 = _tmp117 * _tmp132;
  const Scalar _tmp134 = _tmp121 * _tmp75;
  const Scalar _tmp135 = _tmp132 * _tmp73;
  const Scalar _tmp136 = _tmp122 * _tmp75;
  const Scalar _tmp137 = _tmp132 * _tmp79;
  const Scalar _tmp138 = _tmp130 * _tmp82;
  const Scalar _tmp139 = _tmp123 * _tmp74;
  const Scalar _tmp140 = _tmp139 * _tmp81;
  const Scalar _tmp141 =
      -_tmp122 * _tmp140 + _tmp138 * _tmp81 + _tmp83 * (-_tmp136 * _tmp58 + _tmp137 * _tmp58);
  const Scalar _tmp142 = _tmp60 * _tmp88;
  const Scalar _tmp143 = -_tmp134 * _tmp84 + _tmp135 * _tmp84 - _tmp141 * _tmp85;
  const Scalar _tmp144 = _tmp139 * _tmp99;
  const Scalar _tmp145 =
      -_tmp122 * _tmp144 + _tmp138 * _tmp99 + _tmp83 * (-_tmp136 * _tmp98 + _tmp137 * _tmp98);
  const Scalar _tmp146 = -_tmp100 * _tmp134 + _tmp100 * _tmp135 - _tmp145 * _tmp85;
  const Scalar _tmp147 = _tmp111 * _tmp114;
  const Scalar _tmp148 = _tmp49 * _tmp73;
  const Scalar _tmp149 = _tmp49 * fh1;
  const Scalar _tmp150 = _tmp103 * _tmp149;
  const Scalar _tmp151 = _tmp149 * _tmp95;
  const Scalar _tmp152 = _tmp106 * _tmp109;
  const Scalar _tmp153 = std::exp(-_tmp102 * _tmp150 + _tmp147 * _tmp148 - _tmp151 * _tmp86 +
                                  _tmp152 * _tmp49 * _tmp85);
  const Scalar _tmp154 = _tmp114 * _tmp125;
  const Scalar _tmp155 = _tmp131 * _tmp152;
  const Scalar _tmp156 = _tmp130 * _tmp155;
  const Scalar _tmp157 = _tmp152 * _tmp75;
  const Scalar _tmp158 = _tmp157 * _tmp49;
  const Scalar _tmp159 = _tmp147 * _tmp49;
  const Scalar _tmp160 = _tmp75 * fh1;
  const Scalar _tmp161 = _tmp160 * _tmp95;
  const Scalar _tmp162 = _tmp103 * _tmp160;
  const Scalar _tmp163 = std::exp(-_tmp100 * _tmp162 - _tmp147 - _tmp157 - _tmp161 * _tmp84);
  const Scalar _tmp164 = _tmp132 * fh1;
  const Scalar _tmp165 = _tmp84 * _tmp95;
  const Scalar _tmp166 = _tmp100 * _tmp103;
  const Scalar _tmp167 = _tmp120 * _tmp75;
  const Scalar _tmp168 = _tmp118 * _tmp64;
  const Scalar _tmp169 = _tmp120 * _tmp29 - _tmp168 + _tmp66;
  const Scalar _tmp170 =
      _tmp120 * _tmp38 - _tmp129 - _tmp169 * _tmp57 -
      _tmp56 * (-_tmp120 * _tmp39 + _tmp168 * _tmp69 - _tmp169 * _tmp53 - _tmp70);
  const Scalar _tmp171 = _tmp131 * _tmp170;
  const Scalar _tmp172 = _tmp117 * _tmp171;
  const Scalar _tmp173 =
      _tmp120 * _tmp33 + _tmp120 * _tmp62 - _tmp168 * _tmp59 - _tmp169 * _tmp77 + _tmp59 * _tmp66;
  const Scalar _tmp174 = _tmp173 * _tmp75;
  const Scalar _tmp175 = _tmp171 * _tmp79;
  const Scalar _tmp176 = _tmp170 * _tmp82;
  const Scalar _tmp177 =
      -_tmp144 * _tmp173 + _tmp176 * _tmp99 + _tmp83 * (-_tmp174 * _tmp98 + _tmp175 * _tmp98);
  const Scalar _tmp178 = _tmp169 * _tmp75;
  const Scalar _tmp179 = _tmp171 * _tmp73;
  const Scalar _tmp180 = -_tmp100 * _tmp178 + _tmp100 * _tmp179 - _tmp177 * _tmp85;
  const Scalar _tmp181 = _tmp124 * _tmp173;
  const Scalar _tmp182 =
      -_tmp140 * _tmp173 + _tmp176 * _tmp81 + _tmp83 * (-_tmp174 * _tmp58 + _tmp175 * _tmp58);
  const Scalar _tmp183 = -_tmp178 * _tmp84 + _tmp179 * _tmp84 - _tmp182 * _tmp85;
  const Scalar _tmp184 = _tmp114 * _tmp181;
  const Scalar _tmp185 = _tmp155 * _tmp170;
  const Scalar _tmp186 = _tmp171 * fh1;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) = -_tmp116 * (-_tmp104 * (_tmp100 * _tmp127 - _tmp100 * _tmp133 - _tmp101 +
                                       _tmp142 * _tmp145 + _tmp146 * _tmp87) -
                           _tmp110 * (_tmp106 * _tmp127 - _tmp106 * _tmp133 - _tmp107 -
                                      _tmp108 * _tmp134 + _tmp108 * _tmp135) -
                           _tmp115 * (_tmp111 * _tmp119 - _tmp112 - _tmp113 * _tmp121 -
                                      _tmp117 * _tmp125 + _tmp125 * _tmp126) -
                           _tmp97 * (_tmp127 * _tmp84 - _tmp133 * _tmp84 + _tmp141 * _tmp142 +
                                     _tmp143 * _tmp87 - _tmp89));
  _res(2, 0) = -_tmp153 * (-_tmp121 * _tmp158 - _tmp121 * _tmp159 + _tmp143 * _tmp151 +
                           _tmp146 * _tmp150 + _tmp148 * _tmp154 + _tmp148 * _tmp156);
  _res(3, 0) = -_tmp163 * (_tmp141 * _tmp161 + _tmp145 * _tmp162 - _tmp154 - _tmp156 -
                           _tmp164 * _tmp165 - _tmp164 * _tmp166);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp116 *
      (-_tmp104 * (_tmp100 * _tmp167 - _tmp100 * _tmp172 + _tmp142 * _tmp177 + _tmp180 * _tmp87) -
       _tmp110 * (_tmp106 * _tmp167 - _tmp106 * _tmp172 - _tmp108 * _tmp178 + _tmp108 * _tmp179) -
       _tmp115 * (_tmp111 * _tmp120 - _tmp113 * _tmp169 - _tmp117 * _tmp181 + _tmp126 * _tmp181) -
       _tmp97 * (_tmp142 * _tmp182 + _tmp167 * _tmp84 - _tmp172 * _tmp84 + _tmp183 * _tmp87));
  _res(2, 1) = -_tmp153 * (_tmp148 * _tmp184 + _tmp148 * _tmp185 + _tmp150 * _tmp180 +
                           _tmp151 * _tmp183 - _tmp158 * _tmp169 - _tmp159 * _tmp169);
  _res(3, 1) = -_tmp163 * (_tmp161 * _tmp182 + _tmp162 * _tmp177 - _tmp165 * _tmp186 -
                           _tmp166 * _tmp186 - _tmp184 - _tmp185);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
