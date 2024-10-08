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
 * Symbolic function: IK_residual_func_cost3_wrt_pb_Nl16
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPbNl16(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 569

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (184)
  const Scalar _tmp0 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp6 = 2 * _tmp0;
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp9 = _tmp2 * _tmp8;
  const Scalar _tmp10 =
      -Scalar(0.010999999999999999) * _tmp7 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = 2 * _tmp2 * _tmp5;
  const Scalar _tmp12 = _tmp0 * _tmp8;
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = _tmp10 - _tmp13;
  const Scalar _tmp15 = _tmp14 + _tmp4;
  const Scalar _tmp16 = _tmp15 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp17 = std::pow(_tmp16, Scalar(2));
  const Scalar _tmp18 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp18 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp20 = -_tmp19;
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp22 = _tmp2 * _tmp6;
  const Scalar _tmp23 = _tmp5 * _tmp8;
  const Scalar _tmp24 =
      -Scalar(0.010999999999999999) * _tmp22 + Scalar(0.010999999999999999) * _tmp23;
  const Scalar _tmp25 = _tmp21 + _tmp24;
  const Scalar _tmp26 = _tmp20 + _tmp25;
  const Scalar _tmp27 = _tmp26 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp28 = std::pow(_tmp27, Scalar(2));
  const Scalar _tmp29 = _tmp17 + _tmp28;
  const Scalar _tmp30 = std::pow(_tmp29, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp31 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp32 =
      -Scalar(0.010999999999999999) * _tmp18 - Scalar(0.010999999999999999) * _tmp3;
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp7 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp34 = _tmp32 - _tmp33;
  const Scalar _tmp35 = _tmp31 + _tmp34;
  const Scalar _tmp36 = _tmp30 * _tmp35;
  const Scalar _tmp37 = _tmp16 * _tmp36;
  const Scalar _tmp38 = -_tmp31;
  const Scalar _tmp39 = _tmp32 + _tmp33 + _tmp38;
  const Scalar _tmp40 = _tmp30 * _tmp39;
  const Scalar _tmp41 = -_tmp4;
  const Scalar _tmp42 = _tmp10 + _tmp13;
  const Scalar _tmp43 = _tmp41 + _tmp42;
  const Scalar _tmp44 = _tmp43 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp45 = Scalar(1.0) / (_tmp44);
  const Scalar _tmp46 = -_tmp21 + _tmp24;
  const Scalar _tmp47 = _tmp19 + _tmp46;
  const Scalar _tmp48 = _tmp47 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp49 = _tmp45 * _tmp48;
  const Scalar _tmp50 = _tmp30 * _tmp49;
  const Scalar _tmp51 = _tmp27 * _tmp30;
  const Scalar _tmp52 = _tmp16 * _tmp50 - _tmp51;
  const Scalar _tmp53 = _tmp14 + _tmp41;
  const Scalar _tmp54 = _tmp53 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp55 = _tmp20 + _tmp46;
  const Scalar _tmp56 = _tmp55 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp57 = std::pow(Scalar(std::pow(_tmp54, Scalar(2)) + std::pow(_tmp56, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp58 = _tmp54 * _tmp57;
  const Scalar _tmp59 = _tmp56 * _tmp57;
  const Scalar _tmp60 = Scalar(1.0) / (_tmp49 * _tmp58 - _tmp59);
  const Scalar _tmp61 = _tmp35 * _tmp58;
  const Scalar _tmp62 = _tmp34 + _tmp38;
  const Scalar _tmp63 = _tmp60 * (-_tmp49 * _tmp61 + _tmp59 * _tmp62);
  const Scalar _tmp64 = Scalar(1.0) * _tmp47;
  const Scalar _tmp65 = Scalar(1.0) * _tmp43;
  const Scalar _tmp66 = (-_tmp53 + _tmp65) / (_tmp55 - _tmp64);
  const Scalar _tmp67 = _tmp60 * (-_tmp58 * _tmp62 + _tmp61);
  const Scalar _tmp68 = -_tmp16 * _tmp40 + _tmp37 - _tmp52 * _tmp67 -
                        _tmp66 * (_tmp27 * _tmp40 - _tmp37 * _tmp49 - _tmp52 * _tmp63);
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp64 * _tmp66 + _tmp65;
  const Scalar _tmp71 = 0;
  const Scalar _tmp72 = _tmp69 * _tmp71;
  const Scalar _tmp73 = _tmp30 * _tmp72;
  const Scalar _tmp74 = _tmp58 * _tmp60;
  const Scalar _tmp75 = _tmp72 * _tmp74;
  const Scalar _tmp76 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp77 =
      std::sqrt(Scalar(std::pow(_tmp44, Scalar(2)) + std::pow(_tmp48, Scalar(2))));
  const Scalar _tmp78 = _tmp45 * _tmp77;
  const Scalar _tmp79 = _tmp76 * _tmp78;
  const Scalar _tmp80 = _tmp26 * _tmp30;
  const Scalar _tmp81 = Scalar(1.0) / (_tmp77);
  const Scalar _tmp82 = _tmp78 * (_tmp43 * _tmp48 * _tmp81 - _tmp44 * _tmp47 * _tmp81);
  const Scalar _tmp83 = _tmp60 * (-_tmp53 * _tmp59 + _tmp55 * _tmp58 + _tmp58 * _tmp82);
  const Scalar _tmp84 = _tmp30 * _tmp82;
  const Scalar _tmp85 = -_tmp15 * _tmp51 + _tmp16 * _tmp80 + _tmp16 * _tmp84 - _tmp52 * _tmp83;
  const Scalar _tmp86 = Scalar(1.0) / (_tmp85);
  const Scalar _tmp87 = Scalar(1.0) * _tmp86;
  const Scalar _tmp88 = _tmp30 * _tmp87;
  const Scalar _tmp89 = _tmp74 * _tmp87;
  const Scalar _tmp90 = _tmp19 + _tmp25;
  const Scalar _tmp91 = _tmp4 + _tmp42;
  const Scalar _tmp92 = _tmp91 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp93 = _tmp90 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp94 = std::pow(Scalar(std::pow(_tmp92, Scalar(2)) + std::pow(_tmp93, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp95 = _tmp92 * _tmp94;
  const Scalar _tmp96 = _tmp93 * _tmp94;
  const Scalar _tmp97 = fh1 * (-_tmp90 * _tmp95 + _tmp91 * _tmp96);
  const Scalar _tmp98 = _tmp78 * _tmp97;
  const Scalar _tmp99 = _tmp35 * _tmp49;
  const Scalar _tmp100 = -_tmp35 + _tmp49 * _tmp67 - _tmp66 * (_tmp49 * _tmp63 + _tmp99);
  const Scalar _tmp101 = _tmp69 * _tmp85;
  const Scalar _tmp102 = -_tmp100 * _tmp101 + _tmp49 * _tmp83 - _tmp82;
  const Scalar _tmp103 = _tmp68 * _tmp86;
  const Scalar _tmp104 = _tmp100 + _tmp102 * _tmp103;
  const Scalar _tmp105 = _tmp30 * _tmp69;
  const Scalar _tmp106 = _tmp104 * _tmp105;
  const Scalar _tmp107 = _tmp52 * _tmp69;
  const Scalar _tmp108 = -_tmp104 * _tmp107 - _tmp49;
  const Scalar _tmp109 = _tmp95 * fh1;
  const Scalar _tmp110 = _tmp109 * _tmp78;
  const Scalar _tmp111 = Scalar(1.0) * _tmp63 * _tmp66 - Scalar(1.0) * _tmp67;
  const Scalar _tmp112 = -_tmp101 * _tmp111 - Scalar(1.0) * _tmp83;
  const Scalar _tmp113 = _tmp103 * _tmp112 + _tmp111;
  const Scalar _tmp114 = _tmp105 * _tmp113;
  const Scalar _tmp115 = -_tmp107 * _tmp113 + Scalar(1.0);
  const Scalar _tmp116 = _tmp96 * fh1;
  const Scalar _tmp117 = _tmp116 * _tmp78;
  const Scalar _tmp118 = std::exp(_tmp110 * (_tmp106 * _tmp16 + _tmp108 * _tmp74 + Scalar(1.0)) +
                                  _tmp117 * (_tmp114 * _tmp16 + _tmp115 * _tmp74) +
                                  _tmp79 * (_tmp16 * _tmp73 - _tmp52 * _tmp75) +
                                  _tmp98 * (_tmp16 * _tmp88 - _tmp52 * _tmp89));
  const Scalar _tmp119 = _tmp16 * _tmp30;
  const Scalar _tmp120 = std::pow(_tmp29, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp121 = _tmp120 * _tmp17;
  const Scalar _tmp122 = _tmp121 * _tmp35;
  const Scalar _tmp123 = _tmp120 * _tmp16 * _tmp27;
  const Scalar _tmp124 = _tmp121 * _tmp49 - _tmp123 - _tmp50;
  const Scalar _tmp125 = _tmp123 * _tmp39;
  const Scalar _tmp126 =
      -_tmp121 * _tmp39 + _tmp122 - _tmp124 * _tmp67 - _tmp36 + _tmp40 -
      _tmp66 * (-_tmp122 * _tmp49 - _tmp124 * _tmp63 + _tmp125 + _tmp36 * _tmp49);
  const Scalar _tmp127 = std::pow(_tmp68, Scalar(-2));
  const Scalar _tmp128 = _tmp126 * _tmp127;
  const Scalar _tmp129 = _tmp104 * _tmp128;
  const Scalar _tmp130 = _tmp104 * _tmp69;
  const Scalar _tmp131 = _tmp126 * _tmp86;
  const Scalar _tmp132 = _tmp100 * _tmp85;
  const Scalar _tmp133 =
      _tmp121 * _tmp26 + _tmp121 * _tmp82 - _tmp123 * _tmp15 - _tmp124 * _tmp83 - _tmp80 - _tmp84;
  const Scalar _tmp134 = _tmp133 * _tmp69;
  const Scalar _tmp135 = std::pow(_tmp85, Scalar(-2));
  const Scalar _tmp136 = _tmp135 * _tmp68;
  const Scalar _tmp137 = _tmp102 * _tmp136;
  const Scalar _tmp138 =
      _tmp102 * _tmp131 + _tmp103 * (-_tmp100 * _tmp134 + _tmp128 * _tmp132) - _tmp133 * _tmp137;
  const Scalar _tmp139 = _tmp105 * _tmp16;
  const Scalar _tmp140 = _tmp128 * _tmp52;
  const Scalar _tmp141 = _tmp104 * _tmp140 - _tmp107 * _tmp138 - _tmp124 * _tmp130;
  const Scalar _tmp142 = Scalar(1.0) * _tmp135;
  const Scalar _tmp143 = _tmp133 * _tmp142;
  const Scalar _tmp144 = _tmp52 * _tmp74;
  const Scalar _tmp145 = _tmp113 * _tmp69;
  const Scalar _tmp146 = _tmp119 * _tmp128;
  const Scalar _tmp147 = _tmp111 * _tmp85;
  const Scalar _tmp148 = _tmp112 * _tmp136;
  const Scalar _tmp149 =
      _tmp103 * (-_tmp111 * _tmp134 + _tmp128 * _tmp147) + _tmp112 * _tmp131 - _tmp133 * _tmp148;
  const Scalar _tmp150 = -_tmp107 * _tmp149 + _tmp113 * _tmp140 - _tmp124 * _tmp145;
  const Scalar _tmp151 = _tmp71 * _tmp74;
  const Scalar _tmp152 = _tmp60 * fh1;
  const Scalar _tmp153 = _tmp152 * _tmp96;
  const Scalar _tmp154 = _tmp152 * _tmp95;
  const Scalar _tmp155 = _tmp72 * _tmp76;
  const Scalar _tmp156 = _tmp155 * _tmp60;
  const Scalar _tmp157 = _tmp87 * _tmp97;
  const Scalar _tmp158 = _tmp157 * _tmp60;
  const Scalar _tmp159 =
      std::exp(-_tmp108 * _tmp154 - _tmp115 * _tmp153 + _tmp156 * _tmp52 + _tmp158 * _tmp52);
  const Scalar _tmp160 = _tmp143 * _tmp97;
  const Scalar _tmp161 = _tmp52 * _tmp60;
  const Scalar _tmp162 = _tmp71 * _tmp76;
  const Scalar _tmp163 = _tmp162 * _tmp60;
  const Scalar _tmp164 = std::exp(-_tmp109 * _tmp130 - _tmp116 * _tmp145 - _tmp155 - _tmp157);
  const Scalar _tmp165 = _tmp113 * _tmp116;
  const Scalar _tmp166 = _tmp69 * fh1;
  const Scalar _tmp167 = _tmp166 * _tmp95;
  const Scalar _tmp168 = _tmp166 * _tmp96;
  const Scalar _tmp169 = _tmp120 * _tmp28;
  const Scalar _tmp170 = _tmp123 * _tmp49 - _tmp169 + _tmp30;
  const Scalar _tmp171 =
      _tmp123 * _tmp35 - _tmp125 - _tmp170 * _tmp67 -
      _tmp66 * (-_tmp123 * _tmp99 + _tmp169 * _tmp39 - _tmp170 * _tmp63 - _tmp40);
  const Scalar _tmp172 = _tmp127 * _tmp171;
  const Scalar _tmp173 = _tmp119 * _tmp172;
  const Scalar _tmp174 = _tmp172 * _tmp52;
  const Scalar _tmp175 =
      _tmp123 * _tmp26 + _tmp123 * _tmp82 - _tmp15 * _tmp169 + _tmp15 * _tmp30 - _tmp170 * _tmp83;
  const Scalar _tmp176 = _tmp171 * _tmp86;
  const Scalar _tmp177 = _tmp175 * _tmp69;
  const Scalar _tmp178 =
      _tmp102 * _tmp176 + _tmp103 * (-_tmp100 * _tmp177 + _tmp132 * _tmp172) - _tmp137 * _tmp175;
  const Scalar _tmp179 = _tmp104 * _tmp174 - _tmp107 * _tmp178 - _tmp130 * _tmp170;
  const Scalar _tmp180 = _tmp142 * _tmp175;
  const Scalar _tmp181 =
      _tmp103 * (-_tmp111 * _tmp177 + _tmp147 * _tmp172) + _tmp112 * _tmp176 - _tmp148 * _tmp175;
  const Scalar _tmp182 = -_tmp107 * _tmp181 + _tmp113 * _tmp174 - _tmp145 * _tmp170;
  const Scalar _tmp183 = _tmp180 * _tmp97;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) = -_tmp118 * (-_tmp110 * (-_tmp106 - _tmp119 * _tmp129 + _tmp121 * _tmp130 +
                                       _tmp138 * _tmp139 + _tmp141 * _tmp74) -
                           _tmp117 * (-_tmp113 * _tmp146 - _tmp114 + _tmp121 * _tmp145 +
                                      _tmp139 * _tmp149 + _tmp150 * _tmp74) -
                           _tmp79 * (_tmp121 * _tmp72 - _tmp124 * _tmp75 + _tmp140 * _tmp151 -
                                     _tmp146 * _tmp71 - _tmp73) -
                           _tmp98 * (-_tmp119 * _tmp143 + _tmp121 * _tmp87 - _tmp124 * _tmp89 +
                                     _tmp143 * _tmp144 - _tmp88));
  _res(2, 0) = -_tmp159 * (-_tmp124 * _tmp156 - _tmp124 * _tmp158 + _tmp140 * _tmp163 +
                           _tmp141 * _tmp154 + _tmp150 * _tmp153 + _tmp160 * _tmp161);
  _res(3, 0) = -_tmp164 * (-_tmp109 * _tmp129 - _tmp128 * _tmp162 - _tmp128 * _tmp165 +
                           _tmp138 * _tmp167 + _tmp149 * _tmp168 - _tmp160);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp118 *
      (-_tmp110 * (-_tmp104 * _tmp173 + _tmp123 * _tmp130 + _tmp139 * _tmp178 + _tmp179 * _tmp74) -
       _tmp117 * (-_tmp113 * _tmp173 + _tmp123 * _tmp145 + _tmp139 * _tmp181 + _tmp182 * _tmp74) -
       _tmp79 * (_tmp123 * _tmp72 + _tmp151 * _tmp174 - _tmp170 * _tmp75 - _tmp173 * _tmp71) -
       _tmp98 * (-_tmp119 * _tmp180 + _tmp123 * _tmp87 + _tmp144 * _tmp180 - _tmp170 * _tmp89));
  _res(2, 1) = -_tmp159 * (_tmp153 * _tmp182 + _tmp154 * _tmp179 - _tmp156 * _tmp170 -
                           _tmp158 * _tmp170 + _tmp161 * _tmp183 + _tmp163 * _tmp174);
  _res(3, 1) = -_tmp164 * (-_tmp104 * _tmp109 * _tmp172 - _tmp162 * _tmp172 - _tmp165 * _tmp172 +
                           _tmp167 * _tmp178 + _tmp168 * _tmp181 - _tmp183);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
