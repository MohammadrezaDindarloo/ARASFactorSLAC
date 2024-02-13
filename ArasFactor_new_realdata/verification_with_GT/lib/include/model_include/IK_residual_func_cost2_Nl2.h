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
 * Symbolic function: IK_residual_func_cost2_Nl2
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     p_init0: Scalar
 *     p_init1: Scalar
 *     p_init2: Scalar
 *     rot_init_x: Scalar
 *     rot_init_y: Scalar
 *     rot_init_z: Scalar
 *     rot_init_w: Scalar
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2Nl2(const Scalar fh1, const Scalar fv1,
                                                   const sym::Rot3<Scalar>& DeltaRot,
                                                   const Scalar p_init0, const Scalar p_init1,
                                                   const Scalar p_init2, const Scalar rot_init_x,
                                                   const Scalar rot_init_y, const Scalar rot_init_z,
                                                   const Scalar rot_init_w, const Scalar epsilon) {
  // Total ops: 532

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (161)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp6 = -_tmp5;
  const Scalar _tmp7 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp8 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                       2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp10 = _tmp1 * _tmp9;
  const Scalar _tmp11 = _tmp10 + _tmp8;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = 2 * _tmp1;
  const Scalar _tmp14 = _tmp13 * _tmp7;
  const Scalar _tmp15 = _tmp3 * _tmp9;
  const Scalar _tmp16 = Scalar(0.20999999999999999) * _tmp14 - Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp17 = _tmp12 - _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp6;
  const Scalar _tmp19 = _tmp18 + p_init0;
  const Scalar _tmp20 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp22 = -_tmp21;
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp24 = _tmp13 * _tmp3;
  const Scalar _tmp25 = _tmp7 * _tmp9;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = -Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = -_tmp23 + _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _tmp29 + p_init1;
  const Scalar _tmp31 = Scalar(3.9500536956656402) *
                            std::pow(Scalar(-Scalar(0.50315118556004401) * _tmp19 - 1), Scalar(2)) +
                        Scalar(69.216682114881593) *
                            std::pow(Scalar(-Scalar(0.12019727204189803) * _tmp30 - 1), Scalar(2));
  const Scalar _tmp32 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp33 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp34 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp35 = -_tmp34;
  const Scalar _tmp36 =
      -Scalar(0.010999999999999999) * _tmp2 - Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp37 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp38 = _tmp36 - _tmp37;
  const Scalar _tmp39 = _tmp35 + _tmp38;
  const Scalar _tmp40 = _tmp12 + _tmp16;
  const Scalar _tmp41 = _tmp40 + _tmp5;
  const Scalar _tmp42 = _tmp41 + p_init0;
  const Scalar _tmp43 = _tmp23 + _tmp27;
  const Scalar _tmp44 = _tmp21 + _tmp43;
  const Scalar _tmp45 = _tmp44 + p_init1;
  const Scalar _tmp46 = Scalar(7.3875128562042027) *
                            std::pow(Scalar(1 - Scalar(0.36791786395571047) * _tmp42), Scalar(2)) +
                        Scalar(22.802596067096832) *
                            std::pow(Scalar(1 - Scalar(0.20941503221602112) * _tmp45), Scalar(2));
  const Scalar _tmp47 = _tmp34 + _tmp38;
  const Scalar _tmp48 = _tmp17 + _tmp5;
  const Scalar _tmp49 = _tmp48 + p_init0;
  const Scalar _tmp50 = _tmp49 + Scalar(-2.5202214700000001);
  const Scalar _tmp51 = _tmp22 + _tmp43;
  const Scalar _tmp52 = _tmp51 + p_init1;
  const Scalar _tmp53 = _tmp52 + Scalar(8.3888750099999996);
  const Scalar _tmp54 = std::pow(Scalar(std::pow(_tmp50, Scalar(2)) + std::pow(_tmp53, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp55 = _tmp53 * _tmp54;
  const Scalar _tmp56 = _tmp42 + Scalar(-2.71799795);
  const Scalar _tmp57 = Scalar(1.0) / (_tmp56);
  const Scalar _tmp58 = _tmp45 + Scalar(-4.7752063900000001);
  const Scalar _tmp59 = _tmp57 * _tmp58;
  const Scalar _tmp60 = _tmp36 + _tmp37;
  const Scalar _tmp61 = _tmp34 + _tmp60;
  const Scalar _tmp62 = _tmp50 * _tmp54;
  const Scalar _tmp63 = _tmp61 * _tmp62;
  const Scalar _tmp64 = _tmp47 * _tmp55 - _tmp59 * _tmp63;
  const Scalar _tmp65 = Scalar(1.0) / (-_tmp55 + _tmp59 * _tmp62);
  const Scalar _tmp66 = Scalar(1.0) * _tmp65;
  const Scalar _tmp67 = _tmp64 * _tmp66;
  const Scalar _tmp68 = -_tmp47 * _tmp62 + _tmp63;
  const Scalar _tmp69 = Scalar(1.0) * _tmp44;
  const Scalar _tmp70 = -_tmp69;
  const Scalar _tmp71 = Scalar(1.0) / (_tmp51 + _tmp70);
  const Scalar _tmp72 = Scalar(1.0) * _tmp41;
  const Scalar _tmp73 = _tmp71 * (-_tmp48 + _tmp72);
  const Scalar _tmp74 = -_tmp66 * _tmp68 + _tmp67 * _tmp73;
  const Scalar _tmp75 = _tmp21 + _tmp28;
  const Scalar _tmp76 = _tmp75 + p_init1;
  const Scalar _tmp77 = _tmp76 + Scalar(-4.8333311099999996);
  const Scalar _tmp78 = _tmp40 + _tmp6;
  const Scalar _tmp79 = _tmp78 + p_init0;
  const Scalar _tmp80 = _tmp79 + Scalar(1.79662371);
  const Scalar _tmp81 = std::pow(Scalar(std::pow(_tmp77, Scalar(2)) + std::pow(_tmp80, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp82 = _tmp77 * _tmp81;
  const Scalar _tmp83 = _tmp80 * _tmp81;
  const Scalar _tmp84 = _tmp59 * _tmp83 - _tmp82;
  const Scalar _tmp85 = _tmp65 * _tmp84;
  const Scalar _tmp86 = _tmp59 * _tmp61;
  const Scalar _tmp87 = _tmp35 + _tmp60;
  const Scalar _tmp88 = -_tmp64 * _tmp85 + _tmp82 * _tmp87 - _tmp83 * _tmp86;
  const Scalar _tmp89 = _tmp61 * _tmp83 - _tmp68 * _tmp85 - _tmp73 * _tmp88 - _tmp83 * _tmp87;
  const Scalar _tmp90 = Scalar(1.0) / (_tmp89);
  const Scalar _tmp91 =
      std::sqrt(Scalar(std::pow(_tmp56, Scalar(2)) + std::pow(_tmp58, Scalar(2))));
  const Scalar _tmp92 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp93 = _tmp57 * _tmp91;
  const Scalar _tmp94 = _tmp93 * (_tmp41 * _tmp58 * _tmp92 - _tmp44 * _tmp56 * _tmp92);
  const Scalar _tmp95 = -_tmp48 * _tmp55 + _tmp51 * _tmp62 + _tmp62 * _tmp94;
  const Scalar _tmp96 = _tmp75 * _tmp83 - _tmp78 * _tmp82 + _tmp83 * _tmp94 - _tmp85 * _tmp95;
  const Scalar _tmp97 = _tmp90 * _tmp96;
  const Scalar _tmp98 = Scalar(1.0) / (_tmp96);
  const Scalar _tmp99 = _tmp89 * _tmp98;
  const Scalar _tmp100 = _tmp99 * (-_tmp66 * _tmp95 - _tmp74 * _tmp97);
  const Scalar _tmp101 = _tmp100 + _tmp74;
  const Scalar _tmp102 = _tmp88 * _tmp90;
  const Scalar _tmp103 = _tmp70 + _tmp75;
  const Scalar _tmp104 = _tmp103 * _tmp73;
  const Scalar _tmp105 = Scalar(1.0) / (-_tmp104 + _tmp72 - _tmp78);
  const Scalar _tmp106 = _tmp103 * _tmp105;
  const Scalar _tmp107 = _tmp100 * _tmp106 - _tmp101 * _tmp102 - _tmp67;
  const Scalar _tmp108 = Scalar(1.0) * _tmp71;
  const Scalar _tmp109 = Scalar(1.0) * _tmp105;
  const Scalar _tmp110 = _tmp30 + Scalar(8.3196563700000006);
  const Scalar _tmp111 = _tmp19 + Scalar(1.9874742000000001);
  const Scalar _tmp112 =
      std::pow(Scalar(std::pow(_tmp110, Scalar(2)) + std::pow(_tmp111, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp113 = _tmp110 * _tmp112;
  const Scalar _tmp114 = _tmp113 * fh1;
  const Scalar _tmp115 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp116 = _tmp69 * _tmp73 + _tmp72;
  const Scalar _tmp117 = 0;
  const Scalar _tmp118 = _tmp105 * _tmp116;
  const Scalar _tmp119 = -_tmp103 * _tmp118 - _tmp117 * _tmp88 + _tmp70;
  const Scalar _tmp120 = _tmp111 * _tmp112;
  const Scalar _tmp121 = _tmp39 * fh1;
  const Scalar _tmp122 = Scalar(3.29616) * _tmp11 + _tmp120 * _tmp121 + _tmp18 * fv1;
  const Scalar _tmp123 = _tmp103 * _tmp71;
  const Scalar _tmp124 = Scalar(1.0) * _tmp98;
  const Scalar _tmp125 = _tmp109 * _tmp99;
  const Scalar _tmp126 = _tmp103 * _tmp125 - _tmp124 * _tmp88;
  const Scalar _tmp127 = fh1 * (_tmp113 * _tmp18 - _tmp120 * _tmp29);
  const Scalar _tmp128 = _tmp59 * _tmp65;
  const Scalar _tmp129 = _tmp128 * _tmp64 + _tmp86;
  const Scalar _tmp130 = _tmp128 * _tmp68 - _tmp129 * _tmp73 - _tmp61;
  const Scalar _tmp131 = _tmp99 * (_tmp128 * _tmp95 - _tmp130 * _tmp97 - _tmp94);
  const Scalar _tmp132 = _tmp130 + _tmp131;
  const Scalar _tmp133 = -_tmp102 * _tmp132 + _tmp106 * _tmp131 + _tmp129;
  const Scalar _tmp134 = _tmp120 * fh1;
  const Scalar _tmp135 = -_tmp113 * _tmp121 - Scalar(3.29616) * _tmp26 - _tmp29 * fv1;
  const Scalar _tmp136 = _tmp104 * _tmp109 + Scalar(1.0);
  const Scalar _tmp137 = _tmp109 * _tmp73;
  const Scalar _tmp138 = _tmp83 * _tmp90;
  const Scalar _tmp139 = _tmp84 * _tmp90;
  const Scalar _tmp140 = -_tmp101 * _tmp139 + Scalar(1.0);
  const Scalar _tmp141 = _tmp62 * _tmp65;
  const Scalar _tmp142 = _tmp66 * _tmp84 * _tmp98;
  const Scalar _tmp143 = -_tmp132 * _tmp139 - _tmp59;
  const Scalar _tmp144 = -_tmp114 * _tmp93 * (_tmp101 * _tmp138 + _tmp140 * _tmp141) -
                         _tmp115 * _tmp93 * (-_tmp117 * _tmp62 * _tmp85 + _tmp117 * _tmp83) -
                         _tmp127 * _tmp93 * (_tmp124 * _tmp83 - _tmp142 * _tmp62) -
                         _tmp134 * _tmp93 * (_tmp132 * _tmp138 + _tmp141 * _tmp143 + Scalar(1.0));
  const Scalar _tmp145 = Scalar(1.0) / (_tmp144);
  const Scalar _tmp146 = std::asinh(
      _tmp145 * (Scalar(1.0) * _tmp114 * (_tmp100 * _tmp109 - _tmp107 * _tmp108) +
                 Scalar(1.0) * _tmp115 * (-_tmp108 * _tmp119 - _tmp109 * _tmp116 + Scalar(1.0)) +
                 Scalar(1.0) * _tmp122 * (_tmp109 * _tmp123 - _tmp109) +
                 Scalar(1.0) * _tmp127 * (-_tmp108 * _tmp126 + _tmp125) +
                 Scalar(1.0) * _tmp134 * (-_tmp108 * _tmp133 + _tmp109 * _tmp131) +
                 Scalar(1.0) * _tmp135 * (-_tmp108 * _tmp136 + _tmp137)));
  const Scalar _tmp147 = Scalar(9.6622558468725703) * _tmp144;
  const Scalar _tmp148 = Scalar(6.351516257848961) *
                             std::pow(Scalar(1 - Scalar(0.39679052492160538) * _tmp49), Scalar(2)) +
                         Scalar(70.3732239334025) *
                             std::pow(Scalar(-Scalar(0.11920549523123722) * _tmp52 - 1), Scalar(2));
  const Scalar _tmp149 = _tmp115 * _tmp117;
  const Scalar _tmp150 = _tmp65 * fh1;
  const Scalar _tmp151 = _tmp113 * _tmp140 * _tmp150 + _tmp120 * _tmp143 * _tmp150 -
                         _tmp127 * _tmp142 - _tmp149 * _tmp85;
  const Scalar _tmp152 = Scalar(1.0) / (_tmp151);
  const Scalar _tmp153 = _tmp109 * _tmp122;
  const Scalar _tmp154 =
      std::asinh(_tmp152 * (_tmp107 * _tmp114 * _tmp71 + _tmp115 * _tmp119 * _tmp71 -
                            _tmp123 * _tmp153 + _tmp126 * _tmp127 * _tmp71 +
                            _tmp133 * _tmp134 * _tmp71 + _tmp135 * _tmp136 * _tmp71));
  const Scalar _tmp155 = Scalar(9.6622558468725703) * _tmp151;
  const Scalar _tmp156 = Scalar(23.361089618893828) *
                             std::pow(Scalar(1 - Scalar(0.20689664689659551) * _tmp76), Scalar(2)) +
                         Scalar(3.2278567553341642) *
                             std::pow(Scalar(-Scalar(0.55659957866191134) * _tmp79 - 1), Scalar(2));
  const Scalar _tmp157 =
      _tmp101 * _tmp114 * _tmp90 + _tmp124 * _tmp127 + _tmp132 * _tmp134 * _tmp90 + _tmp149;
  const Scalar _tmp158 = Scalar(1.0) / (_tmp157);
  const Scalar _tmp159 =
      std::asinh(_tmp158 * (-_tmp100 * _tmp105 * _tmp114 - _tmp105 * _tmp131 * _tmp134 +
                            _tmp115 * _tmp118 - _tmp125 * _tmp127 - _tmp135 * _tmp137 + _tmp153));
  const Scalar _tmp160 = Scalar(9.6622558468725703) * _tmp157;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp33 * (-std::sinh(Scalar(1.0) * _tmp32) -
                std::sinh(Scalar(0.1034955) * _tmp0 * (-std::sqrt(_tmp31) - _tmp32 * _tmp33))) -
      Scalar(8.4718465799999993) *
          std::sqrt(Scalar(Scalar(0.013932974275675287) * _tmp31 +
                           std::pow(Scalar(-Scalar(0.11803802046660766) * _tmp39 -
                                           Scalar(0.11803802046660766) * p_init2 + 1),
                                    Scalar(2))));
  _res(1, 0) =
      _tmp147 *
          (-std::sinh(Scalar(1.0) * _tmp146) -
           std::sinh(Scalar(0.1034955) * _tmp145 * (-_tmp146 * _tmp147 - std::sqrt(_tmp46)))) -
      Scalar(8.36416322) *
          std::sqrt(Scalar(Scalar(0.014294040284261563) * _tmp46 +
                           std::pow(Scalar(-Scalar(0.1195576860108189) * _tmp61 -
                                           Scalar(0.1195576860108189) * p_init2 + 1),
                                    Scalar(2))));
  _res(2, 0) =
      _tmp155 *
          (-std::sinh(Scalar(1.0) * _tmp154) -
           std::sinh(Scalar(0.1034955) * _tmp152 * (-std::sqrt(_tmp148) - _tmp154 * _tmp155))) -
      Scalar(8.4693136199999994) *
          std::sqrt(Scalar(Scalar(0.013941309530580858) * _tmp148 +
                           std::pow(Scalar(-Scalar(0.11807332268798426) * _tmp47 -
                                           Scalar(0.11807332268798426) * p_init2 + 1),
                                    Scalar(2))));
  _res(3, 0) =
      _tmp160 *
          (-std::sinh(Scalar(1.0) * _tmp159) -
           std::sinh(Scalar(0.1034955) * _tmp158 * (-std::sqrt(_tmp156) - _tmp159 * _tmp160))) -
      Scalar(8.3700199099999999) *
          std::sqrt(Scalar(Scalar(0.01427404356387209) * _tmp156 +
                           std::pow(Scalar(-Scalar(0.11947402882581673) * _tmp87 -
                                           Scalar(0.11947402882581673) * p_init2 + 1),
                                    Scalar(2))));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
