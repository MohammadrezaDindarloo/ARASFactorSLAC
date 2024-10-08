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
 * Symbolic function: IK_residual_func_cost3_wrt_pa_Nl23
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPaNl23(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 564

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (181)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = 2 * _tmp5;
  const Scalar _tmp7 = _tmp4 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = -_tmp8;
  const Scalar _tmp10 = _tmp2 * _tmp5;
  const Scalar _tmp11 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp12 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp14 = 1 - 2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp9;
  const Scalar _tmp18 = _tmp17 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp19 = std::pow(_tmp18, Scalar(2));
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp21 = -_tmp20;
  const Scalar _tmp22 = _tmp0 * _tmp6;
  const Scalar _tmp23 = _tmp2 * _tmp4;
  const Scalar _tmp24 =
      -Scalar(0.010999999999999999) * _tmp22 + Scalar(0.010999999999999999) * _tmp23;
  const Scalar _tmp25 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp27 = _tmp24 - _tmp26;
  const Scalar _tmp28 = _tmp21 + _tmp27;
  const Scalar _tmp29 = _tmp28 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp30 = std::pow(_tmp29, Scalar(2));
  const Scalar _tmp31 = _tmp19 + _tmp30;
  const Scalar _tmp32 = std::pow(_tmp31, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp10 - Scalar(0.20999999999999999) * _tmp11;
  const Scalar _tmp34 = -Scalar(0.010999999999999999) * _tmp13 -
                        Scalar(0.010999999999999999) * _tmp25 + Scalar(-0.010999999999999999);
  const Scalar _tmp35 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp36 = _tmp34 - _tmp35;
  const Scalar _tmp37 = -_tmp33 + _tmp36;
  const Scalar _tmp38 = _tmp32 * _tmp37;
  const Scalar _tmp39 = _tmp33 + _tmp34 + _tmp35;
  const Scalar _tmp40 = _tmp32 * _tmp39;
  const Scalar _tmp41 = _tmp18 * _tmp40;
  const Scalar _tmp42 = _tmp29 * _tmp32;
  const Scalar _tmp43 = _tmp12 + _tmp15;
  const Scalar _tmp44 = _tmp43 + _tmp8;
  const Scalar _tmp45 = _tmp44 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp46 = Scalar(1.0) / (_tmp45);
  const Scalar _tmp47 = _tmp24 + _tmp26;
  const Scalar _tmp48 = _tmp20 + _tmp47;
  const Scalar _tmp49 = _tmp48 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp50 = _tmp46 * _tmp49;
  const Scalar _tmp51 = _tmp32 * _tmp50;
  const Scalar _tmp52 = _tmp18 * _tmp51 - _tmp42;
  const Scalar _tmp53 = _tmp33 + _tmp36;
  const Scalar _tmp54 = _tmp43 + _tmp9;
  const Scalar _tmp55 = _tmp54 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp56 = _tmp20 + _tmp27;
  const Scalar _tmp57 = _tmp56 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp58 = std::pow(Scalar(std::pow(_tmp55, Scalar(2)) + std::pow(_tmp57, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp59 = _tmp57 * _tmp58;
  const Scalar _tmp60 = _tmp55 * _tmp58;
  const Scalar _tmp61 = _tmp39 * _tmp50;
  const Scalar _tmp62 = Scalar(1.0) / (_tmp50 * _tmp60 - _tmp59);
  const Scalar _tmp63 = _tmp62 * (_tmp53 * _tmp59 - _tmp60 * _tmp61);
  const Scalar _tmp64 = Scalar(1.0) * _tmp48;
  const Scalar _tmp65 = Scalar(1.0) * _tmp44;
  const Scalar _tmp66 = (-_tmp54 + _tmp65) / (_tmp56 - _tmp64);
  const Scalar _tmp67 = _tmp62 * (_tmp39 * _tmp60 - _tmp53 * _tmp60);
  const Scalar _tmp68 = -_tmp18 * _tmp38 + _tmp41 - _tmp52 * _tmp67 -
                        _tmp66 * (_tmp37 * _tmp42 - _tmp41 * _tmp50 - _tmp52 * _tmp63);
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp64 * _tmp66 + _tmp65;
  const Scalar _tmp71 = 0;
  const Scalar _tmp72 = _tmp69 * _tmp71;
  const Scalar _tmp73 = _tmp60 * _tmp62;
  const Scalar _tmp74 = _tmp52 * _tmp73;
  const Scalar _tmp75 = _tmp32 * _tmp72;
  const Scalar _tmp76 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp77 =
      std::sqrt(Scalar(std::pow(_tmp45, Scalar(2)) + std::pow(_tmp49, Scalar(2))));
  const Scalar _tmp78 = _tmp46 * _tmp77;
  const Scalar _tmp79 = _tmp76 * _tmp78;
  const Scalar _tmp80 = Scalar(1.0) / (_tmp77);
  const Scalar _tmp81 = _tmp78 * (_tmp44 * _tmp49 * _tmp80 - _tmp45 * _tmp48 * _tmp80);
  const Scalar _tmp82 = _tmp62 * (-_tmp54 * _tmp59 + _tmp56 * _tmp60 + _tmp60 * _tmp81);
  const Scalar _tmp83 = _tmp18 * _tmp32;
  const Scalar _tmp84 = _tmp28 * _tmp32;
  const Scalar _tmp85 = _tmp17 * _tmp32;
  const Scalar _tmp86 = _tmp18 * _tmp84 - _tmp29 * _tmp85 - _tmp52 * _tmp82 + _tmp81 * _tmp83;
  const Scalar _tmp87 = Scalar(1.0) * _tmp63 * _tmp66 - Scalar(1.0) * _tmp67;
  const Scalar _tmp88 = _tmp69 * _tmp87;
  const Scalar _tmp89 = -Scalar(1.0) * _tmp82 - _tmp86 * _tmp88;
  const Scalar _tmp90 = Scalar(1.0) / (_tmp86);
  const Scalar _tmp91 = _tmp68 * _tmp90;
  const Scalar _tmp92 = _tmp87 + _tmp89 * _tmp91;
  const Scalar _tmp93 = _tmp69 * _tmp92;
  const Scalar _tmp94 = -_tmp52 * _tmp93 + Scalar(1.0);
  const Scalar _tmp95 = _tmp32 * _tmp93;
  const Scalar _tmp96 = _tmp16 + _tmp8;
  const Scalar _tmp97 = _tmp96 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp98 = _tmp21 + _tmp47;
  const Scalar _tmp99 = _tmp98 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp100 = std::pow(Scalar(std::pow(_tmp97, Scalar(2)) + std::pow(_tmp99, Scalar(2))),
                                  Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp101 = _tmp100 * _tmp99;
  const Scalar _tmp102 = _tmp101 * fh1;
  const Scalar _tmp103 = _tmp102 * _tmp78;
  const Scalar _tmp104 = Scalar(1.0) * _tmp90;
  const Scalar _tmp105 = _tmp104 * _tmp32;
  const Scalar _tmp106 = _tmp100 * _tmp97;
  const Scalar _tmp107 = fh1 * (_tmp101 * _tmp96 - _tmp106 * _tmp98);
  const Scalar _tmp108 = _tmp107 * _tmp78;
  const Scalar _tmp109 = -_tmp39 + _tmp50 * _tmp67 - _tmp66 * (_tmp50 * _tmp63 + _tmp61);
  const Scalar _tmp110 = _tmp109 * _tmp69;
  const Scalar _tmp111 = -_tmp110 * _tmp86 + _tmp50 * _tmp82 - _tmp81;
  const Scalar _tmp112 = _tmp109 + _tmp111 * _tmp91;
  const Scalar _tmp113 = _tmp52 * _tmp69;
  const Scalar _tmp114 = -_tmp112 * _tmp113 - _tmp50;
  const Scalar _tmp115 = _tmp112 * _tmp69;
  const Scalar _tmp116 = _tmp115 * _tmp32;
  const Scalar _tmp117 = _tmp106 * fh1;
  const Scalar _tmp118 = _tmp117 * _tmp78;
  const Scalar _tmp119 = std::exp(_tmp103 * (_tmp18 * _tmp95 + _tmp73 * _tmp94) +
                                  _tmp108 * (-_tmp104 * _tmp74 + _tmp105 * _tmp18) +
                                  _tmp118 * (_tmp114 * _tmp73 + _tmp116 * _tmp18 + Scalar(1.0)) +
                                  _tmp79 * (_tmp18 * _tmp75 - _tmp72 * _tmp74));
  const Scalar _tmp120 = std::pow(_tmp31, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp121 = _tmp120 * _tmp19;
  const Scalar _tmp122 = _tmp120 * _tmp18 * _tmp29;
  const Scalar _tmp123 = _tmp121 * _tmp50 - _tmp122 - _tmp51;
  const Scalar _tmp124 = _tmp121 * _tmp28 + _tmp121 * _tmp81 - _tmp122 * _tmp17 - _tmp123 * _tmp82 -
                         _tmp32 * _tmp81 - _tmp84;
  const Scalar _tmp125 = std::pow(_tmp86, Scalar(-2));
  const Scalar _tmp126 = Scalar(1.0) * _tmp125;
  const Scalar _tmp127 = _tmp124 * _tmp126;
  const Scalar _tmp128 = _tmp121 * _tmp39;
  const Scalar _tmp129 = _tmp122 * _tmp37;
  const Scalar _tmp130 =
      -_tmp121 * _tmp37 - _tmp123 * _tmp67 + _tmp128 + _tmp38 - _tmp40 -
      _tmp66 * (-_tmp123 * _tmp63 - _tmp128 * _tmp50 + _tmp129 + _tmp40 * _tmp50);
  const Scalar _tmp131 = std::pow(_tmp68, Scalar(-2));
  const Scalar _tmp132 = _tmp130 * _tmp131;
  const Scalar _tmp133 = _tmp132 * _tmp83;
  const Scalar _tmp134 = _tmp112 * _tmp52;
  const Scalar _tmp135 = _tmp132 * _tmp86;
  const Scalar _tmp136 = _tmp130 * _tmp90;
  const Scalar _tmp137 = _tmp125 * _tmp68;
  const Scalar _tmp138 = _tmp111 * _tmp137;
  const Scalar _tmp139 =
      _tmp111 * _tmp136 - _tmp124 * _tmp138 + _tmp91 * (_tmp109 * _tmp135 - _tmp110 * _tmp124);
  const Scalar _tmp140 = _tmp123 * _tmp69;
  const Scalar _tmp141 = -_tmp112 * _tmp140 - _tmp113 * _tmp139 + _tmp132 * _tmp134;
  const Scalar _tmp142 = _tmp69 * _tmp83;
  const Scalar _tmp143 = _tmp71 * _tmp74;
  const Scalar _tmp144 = _tmp137 * _tmp89;
  const Scalar _tmp145 =
      -_tmp124 * _tmp144 + _tmp136 * _tmp89 + _tmp91 * (-_tmp124 * _tmp88 + _tmp135 * _tmp87);
  const Scalar _tmp146 = _tmp52 * _tmp92;
  const Scalar _tmp147 = -_tmp113 * _tmp145 + _tmp132 * _tmp146 - _tmp140 * _tmp92;
  const Scalar _tmp148 = _tmp62 * fh1;
  const Scalar _tmp149 = _tmp101 * _tmp148;
  const Scalar _tmp150 = _tmp106 * _tmp148;
  const Scalar _tmp151 = _tmp72 * _tmp76;
  const Scalar _tmp152 = _tmp151 * _tmp62;
  const Scalar _tmp153 = _tmp104 * _tmp107;
  const Scalar _tmp154 = _tmp52 * _tmp62;
  const Scalar _tmp155 =
      std::exp(-_tmp114 * _tmp150 - _tmp149 * _tmp94 + _tmp152 * _tmp52 + _tmp153 * _tmp154);
  const Scalar _tmp156 = _tmp71 * _tmp76;
  const Scalar _tmp157 = _tmp153 * _tmp62;
  const Scalar _tmp158 = _tmp132 * _tmp156;
  const Scalar _tmp159 = _tmp107 * _tmp127;
  const Scalar _tmp160 = std::exp(-_tmp102 * _tmp93 - _tmp115 * _tmp117 - _tmp151 - _tmp153);
  const Scalar _tmp161 = _tmp102 * _tmp92;
  const Scalar _tmp162 = _tmp102 * _tmp69;
  const Scalar _tmp163 = _tmp112 * _tmp117;
  const Scalar _tmp164 = _tmp120 * _tmp30;
  const Scalar _tmp165 = _tmp122 * _tmp50 - _tmp164 + _tmp32;
  const Scalar _tmp166 =
      _tmp122 * _tmp28 + _tmp122 * _tmp81 - _tmp164 * _tmp17 - _tmp165 * _tmp82 + _tmp85;
  const Scalar _tmp167 = _tmp126 * _tmp166;
  const Scalar _tmp168 = _tmp165 * _tmp73;
  const Scalar _tmp169 =
      _tmp122 * _tmp39 - _tmp129 - _tmp165 * _tmp67 -
      _tmp66 * (-_tmp122 * _tmp61 + _tmp164 * _tmp37 - _tmp165 * _tmp63 - _tmp38);
  const Scalar _tmp170 = _tmp131 * _tmp169;
  const Scalar _tmp171 = _tmp170 * _tmp83;
  const Scalar _tmp172 = _tmp170 * _tmp86;
  const Scalar _tmp173 = _tmp169 * _tmp90;
  const Scalar _tmp174 =
      _tmp111 * _tmp173 - _tmp138 * _tmp166 + _tmp91 * (_tmp109 * _tmp172 - _tmp110 * _tmp166);
  const Scalar _tmp175 = _tmp174 * _tmp69;
  const Scalar _tmp176 = -_tmp113 * _tmp174 - _tmp115 * _tmp165 + _tmp134 * _tmp170;
  const Scalar _tmp177 =
      -_tmp144 * _tmp166 + _tmp173 * _tmp89 + _tmp91 * (-_tmp166 * _tmp88 + _tmp172 * _tmp87);
  const Scalar _tmp178 = -_tmp113 * _tmp177 + _tmp146 * _tmp170 - _tmp165 * _tmp93;
  const Scalar _tmp179 = _tmp107 * _tmp167;
  const Scalar _tmp180 = _tmp156 * _tmp170;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) = -_tmp119 * (-_tmp103 * (_tmp121 * _tmp93 - _tmp133 * _tmp92 + _tmp142 * _tmp145 +
                                       _tmp147 * _tmp73 - _tmp95) -
                           _tmp108 * (_tmp104 * _tmp121 - _tmp104 * _tmp123 * _tmp73 - _tmp105 +
                                      _tmp127 * _tmp74 - _tmp127 * _tmp83) -
                           _tmp118 * (-_tmp112 * _tmp133 + _tmp115 * _tmp121 - _tmp116 +
                                      _tmp139 * _tmp142 + _tmp141 * _tmp73) -
                           _tmp79 * (_tmp121 * _tmp72 + _tmp132 * _tmp143 - _tmp133 * _tmp71 -
                                     _tmp140 * _tmp71 * _tmp73 - _tmp75));
  _res(2, 0) = -_tmp155 * (-_tmp123 * _tmp157 - _tmp140 * _tmp156 * _tmp62 + _tmp141 * _tmp150 +
                           _tmp147 * _tmp149 + _tmp154 * _tmp158 + _tmp154 * _tmp159);
  _res(3, 0) = -_tmp160 * (_tmp117 * _tmp139 * _tmp69 - _tmp132 * _tmp161 - _tmp132 * _tmp163 +
                           _tmp145 * _tmp162 - _tmp158 - _tmp159);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp119 *
      (-_tmp103 * (_tmp122 * _tmp93 + _tmp142 * _tmp177 - _tmp171 * _tmp92 + _tmp178 * _tmp73) -
       _tmp108 * (_tmp104 * _tmp122 - _tmp104 * _tmp168 + _tmp167 * _tmp74 - _tmp167 * _tmp83) -
       _tmp118 * (-_tmp112 * _tmp171 + _tmp115 * _tmp122 + _tmp175 * _tmp83 + _tmp176 * _tmp73) -
       _tmp79 * (_tmp122 * _tmp72 + _tmp143 * _tmp170 - _tmp168 * _tmp72 - _tmp171 * _tmp71));
  _res(2, 1) = -_tmp155 * (_tmp149 * _tmp178 + _tmp150 * _tmp176 - _tmp152 * _tmp165 +
                           _tmp154 * _tmp179 + _tmp154 * _tmp180 - _tmp157 * _tmp165);
  _res(3, 1) = -_tmp160 * (_tmp117 * _tmp175 - _tmp161 * _tmp170 + _tmp162 * _tmp177 -
                           _tmp163 * _tmp170 - _tmp179 - _tmp180);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
