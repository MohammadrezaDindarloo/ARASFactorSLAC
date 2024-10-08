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
 * Symbolic function: IK_residual_func_cost3_wrt_pd_Nl8
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPdNl8(
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
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = _tmp4 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp9 = _tmp1 * _tmp5;
  const Scalar _tmp10 =
      -Scalar(0.010999999999999999) * _tmp8 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp12 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = _tmp10 - _tmp13;
  const Scalar _tmp15 = _tmp14 + _tmp7;
  const Scalar _tmp16 = _tmp15 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp17 = std::pow(_tmp16, Scalar(2));
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp19 = -_tmp18;
  const Scalar _tmp20 = _tmp2 * _tmp4;
  const Scalar _tmp21 = _tmp0 * _tmp5;
  const Scalar _tmp22 =
      -Scalar(0.010999999999999999) * _tmp20 + Scalar(0.010999999999999999) * _tmp21;
  const Scalar _tmp23 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp11 +
                        Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999);
  const Scalar _tmp25 = _tmp22 + _tmp24;
  const Scalar _tmp26 = _tmp19 + _tmp25;
  const Scalar _tmp27 = _tmp26 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp28 = std::pow(_tmp27, Scalar(2));
  const Scalar _tmp29 = _tmp17 + _tmp28;
  const Scalar _tmp30 = std::pow(_tmp29, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp31 = _tmp10 + _tmp13;
  const Scalar _tmp32 = _tmp31 + _tmp7;
  const Scalar _tmp33 = _tmp32 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp34 = _tmp18 + _tmp25;
  const Scalar _tmp35 = _tmp34 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp36 =
      std::sqrt(Scalar(std::pow(_tmp33, Scalar(2)) + std::pow(_tmp35, Scalar(2))));
  const Scalar _tmp37 = Scalar(1.0) / (_tmp36);
  const Scalar _tmp38 = Scalar(1.0) / (_tmp33);
  const Scalar _tmp39 = _tmp36 * _tmp38;
  const Scalar _tmp40 = _tmp39 * (_tmp32 * _tmp35 * _tmp37 - _tmp33 * _tmp34 * _tmp37);
  const Scalar _tmp41 = _tmp16 * _tmp30;
  const Scalar _tmp42 = _tmp15 * _tmp30;
  const Scalar _tmp43 = _tmp26 * _tmp30;
  const Scalar _tmp44 = -_tmp7;
  const Scalar _tmp45 = _tmp14 + _tmp44;
  const Scalar _tmp46 = _tmp45 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp47 = _tmp22 - _tmp24;
  const Scalar _tmp48 = _tmp19 + _tmp47;
  const Scalar _tmp49 = _tmp48 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp50 = std::pow(Scalar(std::pow(_tmp46, Scalar(2)) + std::pow(_tmp49, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp51 = _tmp46 * _tmp50;
  const Scalar _tmp52 = _tmp49 * _tmp50;
  const Scalar _tmp53 = _tmp40 * _tmp51 - _tmp45 * _tmp52 + _tmp48 * _tmp51;
  const Scalar _tmp54 = _tmp35 * _tmp38;
  const Scalar _tmp55 = Scalar(1.0) / (_tmp51 * _tmp54 - _tmp52);
  const Scalar _tmp56 = _tmp30 * _tmp54;
  const Scalar _tmp57 = _tmp16 * _tmp56 - _tmp27 * _tmp30;
  const Scalar _tmp58 = _tmp55 * _tmp57;
  const Scalar _tmp59 = _tmp16 * _tmp43 - _tmp27 * _tmp42 + _tmp40 * _tmp41 - _tmp53 * _tmp58;
  const Scalar _tmp60 = Scalar(1.0) / (_tmp59);
  const Scalar _tmp61 = Scalar(1.0) * _tmp60;
  const Scalar _tmp62 = _tmp30 * _tmp61;
  const Scalar _tmp63 = _tmp51 * _tmp61;
  const Scalar _tmp64 = _tmp31 + _tmp44;
  const Scalar _tmp65 = _tmp64 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp66 = _tmp18 + _tmp47;
  const Scalar _tmp67 = _tmp66 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp68 = std::pow(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp67, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp69 = _tmp67 * _tmp68;
  const Scalar _tmp70 = _tmp65 * _tmp68;
  const Scalar _tmp71 = fh1 * (_tmp64 * _tmp69 - _tmp66 * _tmp70);
  const Scalar _tmp72 = _tmp39 * _tmp71;
  const Scalar _tmp73 = Scalar(0.20999999999999999) * _tmp8 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp74 =
      -Scalar(0.010999999999999999) * _tmp12 - Scalar(0.010999999999999999) * _tmp23;
  const Scalar _tmp75 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp76 = _tmp74 + _tmp75;
  const Scalar _tmp77 = _tmp73 + _tmp76;
  const Scalar _tmp78 = _tmp30 * _tmp77;
  const Scalar _tmp79 = _tmp16 * _tmp78;
  const Scalar _tmp80 = -_tmp73;
  const Scalar _tmp81 = _tmp76 + _tmp80;
  const Scalar _tmp82 = _tmp30 * _tmp81;
  const Scalar _tmp83 = _tmp54 * _tmp77;
  const Scalar _tmp84 = _tmp74 - _tmp75 + _tmp80;
  const Scalar _tmp85 = -_tmp51 * _tmp83 + _tmp52 * _tmp84;
  const Scalar _tmp86 = Scalar(1.0) * _tmp32;
  const Scalar _tmp87 = Scalar(1.0) * _tmp34;
  const Scalar _tmp88 = (-_tmp45 + _tmp86) / (_tmp48 - _tmp87);
  const Scalar _tmp89 = _tmp51 * _tmp77 - _tmp51 * _tmp84;
  const Scalar _tmp90 = -_tmp16 * _tmp82 - _tmp58 * _tmp89 + _tmp79 -
                        _tmp88 * (_tmp27 * _tmp82 - _tmp54 * _tmp79 - _tmp58 * _tmp85);
  const Scalar _tmp91 = Scalar(1.0) / (_tmp90);
  const Scalar _tmp92 = Scalar(1.0) * _tmp55;
  const Scalar _tmp93 = _tmp85 * _tmp88 * _tmp92 - _tmp89 * _tmp92;
  const Scalar _tmp94 = _tmp91 * _tmp93;
  const Scalar _tmp95 = -_tmp53 * _tmp92 - _tmp59 * _tmp94;
  const Scalar _tmp96 = _tmp60 * _tmp95;
  const Scalar _tmp97 = _tmp90 * _tmp96 + _tmp93;
  const Scalar _tmp98 = _tmp91 * _tmp97;
  const Scalar _tmp99 = -_tmp57 * _tmp98 + Scalar(1.0);
  const Scalar _tmp100 = _tmp51 * _tmp55;
  const Scalar _tmp101 = _tmp30 * _tmp98;
  const Scalar _tmp102 = _tmp69 * fh1;
  const Scalar _tmp103 = _tmp102 * _tmp39;
  const Scalar _tmp104 = _tmp86 + _tmp87 * _tmp88;
  const Scalar _tmp105 = 0;
  const Scalar _tmp106 = _tmp105 * _tmp91;
  const Scalar _tmp107 = _tmp106 * _tmp51;
  const Scalar _tmp108 = _tmp106 * _tmp30;
  const Scalar _tmp109 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp110 = _tmp109 * _tmp39;
  const Scalar _tmp111 = _tmp54 * _tmp55;
  const Scalar _tmp112 = _tmp111 * _tmp89 - _tmp77 - _tmp88 * (_tmp111 * _tmp85 + _tmp83);
  const Scalar _tmp113 = _tmp112 * _tmp91;
  const Scalar _tmp114 = _tmp111 * _tmp53 - _tmp113 * _tmp59 - _tmp40;
  const Scalar _tmp115 = _tmp114 * _tmp60;
  const Scalar _tmp116 = _tmp112 + _tmp115 * _tmp90;
  const Scalar _tmp117 = _tmp116 * _tmp91;
  const Scalar _tmp118 = -_tmp117 * _tmp57 - _tmp54;
  const Scalar _tmp119 = _tmp117 * _tmp30;
  const Scalar _tmp120 = _tmp70 * fh1;
  const Scalar _tmp121 = _tmp120 * _tmp39;
  const Scalar _tmp122 = std::exp(_tmp103 * (_tmp100 * _tmp99 + _tmp101 * _tmp16) +
                                  _tmp110 * (-_tmp107 * _tmp58 + _tmp108 * _tmp16) +
                                  _tmp121 * (_tmp100 * _tmp118 + _tmp119 * _tmp16 + Scalar(1.0)) +
                                  _tmp72 * (_tmp16 * _tmp62 - _tmp58 * _tmp63));
  const Scalar _tmp123 = std::pow(_tmp29, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp124 = _tmp123 * _tmp17;
  const Scalar _tmp125 = _tmp123 * _tmp16 * _tmp27;
  const Scalar _tmp126 = _tmp124 * _tmp54 - _tmp125 - _tmp56;
  const Scalar _tmp127 = _tmp126 * _tmp55;
  const Scalar _tmp128 = _tmp124 * _tmp26 + _tmp124 * _tmp40 - _tmp125 * _tmp15 - _tmp127 * _tmp53 -
                         _tmp30 * _tmp40 - _tmp43;
  const Scalar _tmp129 = std::pow(_tmp59, Scalar(-2));
  const Scalar _tmp130 = Scalar(1.0) * _tmp129;
  const Scalar _tmp131 = _tmp128 * _tmp130;
  const Scalar _tmp132 = _tmp51 * _tmp58;
  const Scalar _tmp133 = _tmp125 * _tmp81;
  const Scalar _tmp134 =
      _tmp124 * _tmp77 - _tmp124 * _tmp81 - _tmp127 * _tmp89 - _tmp78 + _tmp82 -
      _tmp88 * (-_tmp124 * _tmp83 - _tmp127 * _tmp85 + _tmp133 + _tmp54 * _tmp78);
  const Scalar _tmp135 = _tmp129 * _tmp90;
  const Scalar _tmp136 = _tmp114 * _tmp135;
  const Scalar _tmp137 = std::pow(_tmp90, Scalar(-2));
  const Scalar _tmp138 = _tmp134 * _tmp137;
  const Scalar _tmp139 = _tmp138 * _tmp59;
  const Scalar _tmp140 = _tmp60 * _tmp90;
  const Scalar _tmp141 =
      _tmp115 * _tmp134 - _tmp128 * _tmp136 + _tmp140 * (_tmp112 * _tmp139 - _tmp113 * _tmp128);
  const Scalar _tmp142 = _tmp41 * _tmp91;
  const Scalar _tmp143 = _tmp138 * _tmp41;
  const Scalar _tmp144 = _tmp116 * _tmp57;
  const Scalar _tmp145 = _tmp57 * _tmp91;
  const Scalar _tmp146 = -_tmp117 * _tmp126 + _tmp138 * _tmp144 - _tmp141 * _tmp145;
  const Scalar _tmp147 = _tmp105 * _tmp132;
  const Scalar _tmp148 = _tmp57 * _tmp97;
  const Scalar _tmp149 = _tmp135 * _tmp95;
  const Scalar _tmp150 =
      -_tmp128 * _tmp149 + _tmp134 * _tmp96 + _tmp140 * (-_tmp128 * _tmp94 + _tmp139 * _tmp93);
  const Scalar _tmp151 = -_tmp126 * _tmp98 + _tmp138 * _tmp148 - _tmp145 * _tmp150;
  const Scalar _tmp152 = _tmp55 * fh1;
  const Scalar _tmp153 = _tmp152 * _tmp70;
  const Scalar _tmp154 = _tmp152 * _tmp69;
  const Scalar _tmp155 = _tmp61 * _tmp71;
  const Scalar _tmp156 = _tmp106 * _tmp109;
  const Scalar _tmp157 =
      std::exp(-_tmp118 * _tmp153 - _tmp154 * _tmp99 + _tmp155 * _tmp58 + _tmp156 * _tmp58);
  const Scalar _tmp158 = _tmp131 * _tmp71;
  const Scalar _tmp159 = _tmp105 * _tmp109;
  const Scalar _tmp160 = _tmp138 * _tmp159;
  const Scalar _tmp161 = std::exp(-_tmp102 * _tmp98 - _tmp117 * _tmp120 - _tmp155 - _tmp156);
  const Scalar _tmp162 = _tmp102 * _tmp97;
  const Scalar _tmp163 = _tmp116 * _tmp120;
  const Scalar _tmp164 = _tmp102 * _tmp91;
  const Scalar _tmp165 = _tmp123 * _tmp28;
  const Scalar _tmp166 = _tmp125 * _tmp54 - _tmp165 + _tmp30;
  const Scalar _tmp167 = _tmp166 * _tmp55;
  const Scalar _tmp168 =
      _tmp125 * _tmp26 + _tmp125 * _tmp40 - _tmp15 * _tmp165 - _tmp167 * _tmp53 + _tmp42;
  const Scalar _tmp169 = _tmp130 * _tmp168;
  const Scalar _tmp170 =
      _tmp125 * _tmp77 - _tmp133 - _tmp167 * _tmp89 -
      _tmp88 * (-_tmp125 * _tmp83 + _tmp165 * _tmp81 - _tmp167 * _tmp85 - _tmp82);
  const Scalar _tmp171 = _tmp137 * _tmp170;
  const Scalar _tmp172 = _tmp171 * _tmp41;
  const Scalar _tmp173 = _tmp171 * _tmp59;
  const Scalar _tmp174 =
      _tmp140 * (-_tmp168 * _tmp94 + _tmp173 * _tmp93) - _tmp149 * _tmp168 + _tmp170 * _tmp96;
  const Scalar _tmp175 = -_tmp145 * _tmp174 + _tmp148 * _tmp171 - _tmp166 * _tmp98;
  const Scalar _tmp176 =
      _tmp115 * _tmp170 - _tmp136 * _tmp168 + _tmp140 * (_tmp112 * _tmp173 - _tmp113 * _tmp168);
  const Scalar _tmp177 = _tmp176 * _tmp91;
  const Scalar _tmp178 = -_tmp117 * _tmp166 + _tmp144 * _tmp171 - _tmp145 * _tmp176;
  const Scalar _tmp179 = _tmp159 * _tmp171;
  const Scalar _tmp180 = _tmp169 * _tmp71;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) = -_tmp122 * (-_tmp103 * (_tmp100 * _tmp151 - _tmp101 + _tmp124 * _tmp98 +
                                       _tmp142 * _tmp150 - _tmp143 * _tmp97) -
                           _tmp110 * (-_tmp105 * _tmp143 + _tmp106 * _tmp124 - _tmp107 * _tmp127 -
                                      _tmp108 + _tmp138 * _tmp147) -
                           _tmp121 * (_tmp100 * _tmp146 - _tmp116 * _tmp143 + _tmp117 * _tmp124 -
                                      _tmp119 + _tmp141 * _tmp142) -
                           _tmp72 * (_tmp124 * _tmp61 - _tmp127 * _tmp63 + _tmp131 * _tmp132 -
                                     _tmp131 * _tmp41 - _tmp62));
  _res(2, 0) = -_tmp157 * (-_tmp127 * _tmp155 - _tmp127 * _tmp156 + _tmp146 * _tmp153 +
                           _tmp151 * _tmp154 + _tmp158 * _tmp58 + _tmp160 * _tmp58);
  _res(3, 0) = -_tmp161 * (_tmp120 * _tmp141 * _tmp91 - _tmp138 * _tmp162 - _tmp138 * _tmp163 +
                           _tmp150 * _tmp164 - _tmp158 - _tmp160);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp122 *
      (-_tmp103 * (_tmp100 * _tmp175 + _tmp125 * _tmp98 + _tmp142 * _tmp174 - _tmp172 * _tmp97) -
       _tmp110 * (-_tmp105 * _tmp172 + _tmp106 * _tmp125 - _tmp107 * _tmp167 + _tmp147 * _tmp171) -
       _tmp121 * (_tmp100 * _tmp178 - _tmp116 * _tmp172 + _tmp117 * _tmp125 + _tmp177 * _tmp41) -
       _tmp72 * (_tmp125 * _tmp61 + _tmp132 * _tmp169 - _tmp167 * _tmp63 - _tmp169 * _tmp41));
  _res(2, 1) = -_tmp157 * (_tmp153 * _tmp178 + _tmp154 * _tmp175 - _tmp155 * _tmp167 -
                           _tmp156 * _tmp167 + _tmp179 * _tmp58 + _tmp180 * _tmp58);
  _res(3, 1) = -_tmp161 * (_tmp120 * _tmp177 - _tmp162 * _tmp171 - _tmp163 * _tmp171 +
                           _tmp164 * _tmp174 - _tmp179 - _tmp180);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
