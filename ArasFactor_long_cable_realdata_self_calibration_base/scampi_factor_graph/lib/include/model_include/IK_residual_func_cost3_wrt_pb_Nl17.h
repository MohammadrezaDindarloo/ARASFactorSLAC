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
 * Symbolic function: IK_residual_func_cost3_wrt_pb_Nl17
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPbNl17(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 665

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (213)
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
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp10 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp11 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp12 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp13 = _tmp1 * _tmp5;
  const Scalar _tmp14 =
      -Scalar(0.010999999999999999) * _tmp12 - Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 + _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp8;
  const Scalar _tmp17 = _tmp16 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp18 = _tmp2 * _tmp4;
  const Scalar _tmp19 = _tmp0 * _tmp5;
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp18 + Scalar(0.20999999999999999) * _tmp19;
  const Scalar _tmp21 = -_tmp20;
  const Scalar _tmp22 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp23 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = _tmp21 + _tmp25;
  const Scalar _tmp27 = -_tmp11 + _tmp14;
  const Scalar _tmp28 = _tmp27 + _tmp8;
  const Scalar _tmp29 = _tmp28 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp30 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp31 = -_tmp30;
  const Scalar _tmp32 =
      -Scalar(0.010999999999999999) * _tmp18 + Scalar(0.010999999999999999) * _tmp19;
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp9 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp34 = _tmp32 - _tmp33;
  const Scalar _tmp35 = _tmp31 + _tmp34;
  const Scalar _tmp36 = _tmp35 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp37 = std::pow(Scalar(std::pow(_tmp29, Scalar(2)) + std::pow(_tmp36, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp38 = _tmp29 * _tmp37;
  const Scalar _tmp39 = _tmp30 + _tmp34;
  const Scalar _tmp40 = _tmp39 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp41 = _tmp21 + _tmp23 + _tmp24;
  const Scalar _tmp42 = std::pow(_tmp17, Scalar(2));
  const Scalar _tmp43 = std::pow(_tmp40, Scalar(2));
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = std::pow(_tmp44, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp46 = _tmp41 * _tmp45;
  const Scalar _tmp47 = _tmp27 + _tmp7;
  const Scalar _tmp48 = _tmp47 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp49 = Scalar(1.0) / (_tmp48);
  const Scalar _tmp50 = _tmp32 + _tmp33;
  const Scalar _tmp51 = _tmp31 + _tmp50;
  const Scalar _tmp52 = _tmp51 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp53 = _tmp49 * _tmp52;
  const Scalar _tmp54 = _tmp20 + _tmp25;
  const Scalar _tmp55 = _tmp45 * _tmp54;
  const Scalar _tmp56 = _tmp17 * _tmp55;
  const Scalar _tmp57 = _tmp40 * _tmp46 - _tmp53 * _tmp56;
  const Scalar _tmp58 = _tmp45 * _tmp53;
  const Scalar _tmp59 = _tmp40 * _tmp45;
  const Scalar _tmp60 = _tmp17 * _tmp58 - _tmp59;
  const Scalar _tmp61 = Scalar(1.0) / (_tmp60);
  const Scalar _tmp62 = _tmp36 * _tmp37;
  const Scalar _tmp63 = _tmp38 * _tmp53 - _tmp62;
  const Scalar _tmp64 = _tmp61 * _tmp63;
  const Scalar _tmp65 = _tmp38 * _tmp54;
  const Scalar _tmp66 = Scalar(1.0) * _tmp47;
  const Scalar _tmp67 = Scalar(1.0) * _tmp51;
  const Scalar _tmp68 = (-_tmp16 + _tmp66) / (_tmp39 - _tmp67);
  const Scalar _tmp69 = -_tmp17 * _tmp46 + _tmp56;
  const Scalar _tmp70 = -_tmp26 * _tmp38 - _tmp64 * _tmp69 + _tmp65 -
                        _tmp68 * (_tmp26 * _tmp62 - _tmp53 * _tmp65 - _tmp57 * _tmp64);
  const Scalar _tmp71 =
      std::sqrt(Scalar(std::pow(_tmp48, Scalar(2)) + std::pow(_tmp52, Scalar(2))));
  const Scalar _tmp72 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp73 = _tmp49 * _tmp71;
  const Scalar _tmp74 = _tmp73 * (_tmp47 * _tmp52 * _tmp72 - _tmp48 * _tmp51 * _tmp72);
  const Scalar _tmp75 = _tmp39 * _tmp45;
  const Scalar _tmp76 = _tmp45 * _tmp74;
  const Scalar _tmp77 = -_tmp16 * _tmp59 + _tmp17 * _tmp75 + _tmp17 * _tmp76;
  const Scalar _tmp78 = -_tmp28 * _tmp62 + _tmp35 * _tmp38 + _tmp38 * _tmp74 - _tmp64 * _tmp77;
  const Scalar _tmp79 = Scalar(1.0) / (_tmp78);
  const Scalar _tmp80 = Scalar(1.0) / (_tmp70);
  const Scalar _tmp81 = Scalar(1.0) * _tmp61;
  const Scalar _tmp82 = _tmp68 * _tmp81;
  const Scalar _tmp83 = _tmp57 * _tmp82 - _tmp69 * _tmp81;
  const Scalar _tmp84 = _tmp80 * _tmp83;
  const Scalar _tmp85 = -_tmp77 * _tmp81 - _tmp78 * _tmp84;
  const Scalar _tmp86 = _tmp79 * _tmp85;
  const Scalar _tmp87 = _tmp70 * _tmp86 + _tmp83;
  const Scalar _tmp88 = _tmp63 * _tmp80;
  const Scalar _tmp89 = -_tmp87 * _tmp88 + Scalar(1.0);
  const Scalar _tmp90 = _tmp45 * _tmp61;
  const Scalar _tmp91 = _tmp89 * _tmp90;
  const Scalar _tmp92 = _tmp38 * _tmp80;
  const Scalar _tmp93 = _tmp30 + _tmp50;
  const Scalar _tmp94 = _tmp93 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp95 = _tmp15 + _tmp7;
  const Scalar _tmp96 = _tmp95 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp97 = std::pow(Scalar(std::pow(_tmp94, Scalar(2)) + std::pow(_tmp96, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp98 = _tmp94 * _tmp97;
  const Scalar _tmp99 = _tmp73 * fh1;
  const Scalar _tmp100 = _tmp98 * _tmp99;
  const Scalar _tmp101 = _tmp17 * _tmp90;
  const Scalar _tmp102 = _tmp66 + _tmp67 * _tmp68;
  const Scalar _tmp103 = 0;
  const Scalar _tmp104 = _tmp103 * _tmp88;
  const Scalar _tmp105 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp106 = _tmp105 * _tmp73;
  const Scalar _tmp107 = _tmp53 * _tmp61;
  const Scalar _tmp108 = _tmp53 * _tmp54;
  const Scalar _tmp109 = _tmp107 * _tmp69 - _tmp54 - _tmp68 * (_tmp107 * _tmp57 + _tmp108);
  const Scalar _tmp110 = _tmp109 * _tmp80;
  const Scalar _tmp111 = _tmp107 * _tmp77 - _tmp110 * _tmp78 - _tmp74;
  const Scalar _tmp112 = _tmp111 * _tmp79;
  const Scalar _tmp113 = _tmp109 + _tmp112 * _tmp70;
  const Scalar _tmp114 = -_tmp113 * _tmp88 - _tmp53;
  const Scalar _tmp115 = _tmp114 * _tmp90;
  const Scalar _tmp116 = _tmp96 * _tmp97;
  const Scalar _tmp117 = _tmp116 * _tmp99;
  const Scalar _tmp118 = _tmp17 * _tmp45;
  const Scalar _tmp119 = Scalar(1.0) * _tmp79;
  const Scalar _tmp120 = _tmp119 * _tmp64;
  const Scalar _tmp121 = fh1 * (-_tmp116 * _tmp93 + _tmp95 * _tmp98);
  const Scalar _tmp122 = _tmp121 * _tmp73;
  const Scalar _tmp123 = std::exp(_tmp100 * (_tmp17 * _tmp91 + _tmp87 * _tmp92) +
                                  _tmp106 * (-_tmp101 * _tmp104 + _tmp103 * _tmp92) +
                                  _tmp117 * (_tmp113 * _tmp92 + _tmp115 * _tmp17 + Scalar(1.0)) +
                                  _tmp122 * (-_tmp118 * _tmp120 + _tmp119 * _tmp38));
  const Scalar _tmp124 = _tmp118 * _tmp63;
  const Scalar _tmp125 = std::pow(_tmp60, Scalar(-2));
  const Scalar _tmp126 = std::pow(_tmp44, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp127 = _tmp126 * _tmp42;
  const Scalar _tmp128 = _tmp126 * _tmp17 * _tmp40;
  const Scalar _tmp129 = _tmp125 * (_tmp127 * _tmp53 - _tmp128 - _tmp58);
  const Scalar _tmp130 = Scalar(1.0) * _tmp129;
  const Scalar _tmp131 = _tmp130 * _tmp79;
  const Scalar _tmp132 = _tmp127 * _tmp39 + _tmp127 * _tmp74 - _tmp128 * _tmp16 - _tmp75 - _tmp76;
  const Scalar _tmp133 = _tmp129 * _tmp63;
  const Scalar _tmp134 = -_tmp132 * _tmp64 + _tmp133 * _tmp77;
  const Scalar _tmp135 = std::pow(_tmp78, Scalar(-2));
  const Scalar _tmp136 = Scalar(1.0) * _tmp135;
  const Scalar _tmp137 = _tmp136 * _tmp38;
  const Scalar _tmp138 = _tmp118 * _tmp64;
  const Scalar _tmp139 = _tmp136 * _tmp138;
  const Scalar _tmp140 = _tmp127 * _tmp61;
  const Scalar _tmp141 = _tmp127 * _tmp54;
  const Scalar _tmp142 = -_tmp127 * _tmp41 + _tmp141 + _tmp46 - _tmp55;
  const Scalar _tmp143 = _tmp128 * _tmp41;
  const Scalar _tmp144 = -_tmp141 * _tmp53 + _tmp143 + _tmp53 * _tmp55;
  const Scalar _tmp145 = _tmp57 * _tmp63;
  const Scalar _tmp146 =
      _tmp133 * _tmp69 - _tmp142 * _tmp64 - _tmp68 * (_tmp129 * _tmp145 - _tmp144 * _tmp64);
  const Scalar _tmp147 = std::pow(_tmp70, Scalar(-2));
  const Scalar _tmp148 = _tmp146 * _tmp147;
  const Scalar _tmp149 = _tmp148 * _tmp78;
  const Scalar _tmp150 = _tmp53 * _tmp57;
  const Scalar _tmp151 = _tmp129 * _tmp53;
  const Scalar _tmp152 =
      _tmp107 * _tmp142 - _tmp151 * _tmp69 - _tmp68 * (_tmp107 * _tmp144 - _tmp129 * _tmp150);
  const Scalar _tmp153 = _tmp78 * _tmp80;
  const Scalar _tmp154 = _tmp70 * _tmp79;
  const Scalar _tmp155 = _tmp135 * _tmp70;
  const Scalar _tmp156 = _tmp111 * _tmp155;
  const Scalar _tmp157 = _tmp112 * _tmp146 - _tmp134 * _tmp156 + _tmp152 +
                         _tmp154 * (_tmp107 * _tmp132 + _tmp109 * _tmp149 - _tmp110 * _tmp134 -
                                    _tmp151 * _tmp77 - _tmp152 * _tmp153);
  const Scalar _tmp158 = _tmp148 * _tmp63;
  const Scalar _tmp159 = _tmp113 * _tmp158 - _tmp157 * _tmp88;
  const Scalar _tmp160 = _tmp118 * _tmp129;
  const Scalar _tmp161 = _tmp113 * _tmp38;
  const Scalar _tmp162 = _tmp155 * _tmp85;
  const Scalar _tmp163 = _tmp57 * _tmp68;
  const Scalar _tmp164 =
      -_tmp130 * _tmp163 + _tmp130 * _tmp69 - _tmp142 * _tmp81 + _tmp144 * _tmp82;
  const Scalar _tmp165 = -_tmp134 * _tmp162 + _tmp146 * _tmp86 +
                         _tmp154 * (_tmp130 * _tmp77 - _tmp132 * _tmp81 - _tmp134 * _tmp84 +
                                    _tmp149 * _tmp83 - _tmp153 * _tmp164) +
                         _tmp164;
  const Scalar _tmp166 = _tmp158 * _tmp87 - _tmp165 * _tmp88;
  const Scalar _tmp167 = _tmp148 * _tmp38;
  const Scalar _tmp168 = _tmp103 * _tmp138;
  const Scalar _tmp169 = _tmp61 * fh1;
  const Scalar _tmp170 = _tmp116 * _tmp169;
  const Scalar _tmp171 = _tmp103 * _tmp105;
  const Scalar _tmp172 = _tmp171 * _tmp88;
  const Scalar _tmp173 = _tmp119 * _tmp121;
  const Scalar _tmp174 = _tmp169 * _tmp98;
  const Scalar _tmp175 =
      std::exp(-_tmp114 * _tmp170 + _tmp172 * _tmp61 + _tmp173 * _tmp64 - _tmp174 * _tmp89);
  const Scalar _tmp176 = _tmp148 * _tmp171;
  const Scalar _tmp177 = _tmp116 * fh1;
  const Scalar _tmp178 = _tmp98 * fh1;
  const Scalar _tmp179 = _tmp121 * _tmp136;
  const Scalar _tmp180 = _tmp134 * _tmp179;
  const Scalar _tmp181 = _tmp80 * fh1;
  const Scalar _tmp182 = _tmp181 * _tmp98;
  const Scalar _tmp183 = _tmp116 * _tmp181;
  const Scalar _tmp184 =
      std::exp(-_tmp113 * _tmp183 - _tmp171 * _tmp80 - _tmp173 - _tmp182 * _tmp87);
  const Scalar _tmp185 = _tmp113 * _tmp116;
  const Scalar _tmp186 = _tmp148 * fh1;
  const Scalar _tmp187 = _tmp87 * _tmp98;
  const Scalar _tmp188 = _tmp128 * _tmp61;
  const Scalar _tmp189 = _tmp126 * _tmp43;
  const Scalar _tmp190 = -_tmp108 * _tmp128 + _tmp189 * _tmp41 - _tmp46;
  const Scalar _tmp191 = _tmp125 * (_tmp128 * _tmp53 - _tmp189 + _tmp45);
  const Scalar _tmp192 = _tmp191 * _tmp69;
  const Scalar _tmp193 = _tmp128 * _tmp54 - _tmp143;
  const Scalar _tmp194 =
      _tmp192 * _tmp63 - _tmp193 * _tmp64 - _tmp68 * (_tmp145 * _tmp191 - _tmp190 * _tmp64);
  const Scalar _tmp195 = _tmp147 * _tmp194;
  const Scalar _tmp196 = _tmp195 * _tmp38;
  const Scalar _tmp197 = _tmp191 * _tmp77;
  const Scalar _tmp198 =
      _tmp61 * (_tmp128 * _tmp39 + _tmp128 * _tmp74 - _tmp16 * _tmp189 + _tmp16 * _tmp45);
  const Scalar _tmp199 = _tmp197 * _tmp63 - _tmp198 * _tmp63;
  const Scalar _tmp200 = -Scalar(1.0) * _tmp163 * _tmp191 + _tmp190 * _tmp82 +
                         Scalar(1.0) * _tmp192 - _tmp193 * _tmp81;
  const Scalar _tmp201 = _tmp195 * _tmp78;
  const Scalar _tmp202 = _tmp154 * (-_tmp153 * _tmp200 + Scalar(1.0) * _tmp197 -
                                    Scalar(1.0) * _tmp198 - _tmp199 * _tmp84 + _tmp201 * _tmp83) -
                         _tmp162 * _tmp199 + _tmp194 * _tmp86 + _tmp200;
  const Scalar _tmp203 = _tmp195 * _tmp63;
  const Scalar _tmp204 = -_tmp202 * _tmp88 + _tmp203 * _tmp87;
  const Scalar _tmp205 = _tmp191 * _tmp89;
  const Scalar _tmp206 =
      _tmp107 * _tmp193 - _tmp192 * _tmp53 - _tmp68 * (_tmp107 * _tmp190 - _tmp150 * _tmp191);
  const Scalar _tmp207 = _tmp112 * _tmp194 +
                         _tmp154 * (_tmp109 * _tmp201 - _tmp110 * _tmp199 - _tmp153 * _tmp206 -
                                    _tmp197 * _tmp53 + _tmp198 * _tmp53) -
                         _tmp156 * _tmp199 + _tmp206;
  const Scalar _tmp208 = _tmp114 * _tmp191;
  const Scalar _tmp209 = _tmp113 * _tmp203 - _tmp207 * _tmp88;
  const Scalar _tmp210 = _tmp171 * _tmp195;
  const Scalar _tmp211 = _tmp179 * _tmp199;
  const Scalar _tmp212 = _tmp195 * fh1;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) = -_tmp123 * (-_tmp100 * (_tmp101 * _tmp166 + _tmp140 * _tmp89 - _tmp160 * _tmp89 +
                                       _tmp165 * _tmp92 - _tmp167 * _tmp87 - _tmp91) -
                           _tmp106 * (-_tmp103 * _tmp167 - _tmp104 * _tmp140 + _tmp104 * _tmp160 +
                                      _tmp104 * _tmp90 + _tmp148 * _tmp168) -
                           _tmp117 * (_tmp101 * _tmp159 + _tmp114 * _tmp140 - _tmp114 * _tmp160 -
                                      _tmp115 - _tmp148 * _tmp161 + _tmp157 * _tmp92) -
                           _tmp122 * (-_tmp120 * _tmp127 + _tmp120 * _tmp45 + _tmp124 * _tmp131 -
                                      _tmp134 * _tmp137 + _tmp134 * _tmp139));
  _res(2, 0) = -_tmp175 * (-_tmp114 * _tmp129 * _tmp177 + _tmp121 * _tmp131 * _tmp63 +
                           _tmp129 * _tmp172 - _tmp129 * _tmp178 * _tmp89 + _tmp159 * _tmp170 +
                           _tmp166 * _tmp174 + _tmp176 * _tmp64 + _tmp180 * _tmp64);
  _res(3, 0) = -_tmp184 * (_tmp157 * _tmp183 + _tmp165 * _tmp182 - _tmp176 - _tmp180 -
                           _tmp185 * _tmp186 - _tmp186 * _tmp187);
  _res(0, 1) = 0;
  _res(1, 1) = -_tmp123 * (-_tmp100 * (_tmp101 * _tmp204 - _tmp118 * _tmp205 + _tmp188 * _tmp89 -
                                       _tmp196 * _tmp87 + _tmp202 * _tmp92) -
                           _tmp106 * (-_tmp103 * _tmp196 + _tmp104 * _tmp118 * _tmp191 -
                                      _tmp104 * _tmp188 + _tmp168 * _tmp195) -
                           _tmp117 * (_tmp101 * _tmp209 + _tmp114 * _tmp188 - _tmp118 * _tmp208 -
                                      _tmp161 * _tmp195 + _tmp207 * _tmp92) -
                           _tmp122 * (_tmp119 * _tmp124 * _tmp191 - _tmp120 * _tmp128 -
                                      _tmp137 * _tmp199 + _tmp139 * _tmp199));
  _res(2, 1) = -_tmp175 * (_tmp170 * _tmp209 + _tmp172 * _tmp191 + _tmp173 * _tmp191 * _tmp63 +
                           _tmp174 * _tmp204 - _tmp177 * _tmp208 - _tmp178 * _tmp205 +
                           _tmp210 * _tmp64 + _tmp211 * _tmp64);
  _res(3, 1) = -_tmp184 * (_tmp182 * _tmp202 + _tmp183 * _tmp207 - _tmp185 * _tmp212 -
                           _tmp187 * _tmp212 - _tmp210 - _tmp211);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
