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
 * Symbolic function: IK_residual_func_cost2_Nl21
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2Nl21(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 533

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (160)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp5 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp6 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp7 = -_tmp6;
  const Scalar _tmp8 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp9 = 2 * _tmp2;
  const Scalar _tmp10 = _tmp8 * _tmp9;
  const Scalar _tmp11 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                        2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp12 = _tmp11 * _tmp4;
  const Scalar _tmp13 = _tmp10 + _tmp12;
  const Scalar _tmp14 = -Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp15 = 2 * _tmp4 * _tmp8;
  const Scalar _tmp16 = _tmp11 * _tmp2;
  const Scalar _tmp17 = Scalar(0.20999999999999999) * _tmp15 - Scalar(0.20999999999999999) * _tmp16;
  const Scalar _tmp18 = _tmp14 + _tmp17;
  const Scalar _tmp19 = _tmp18 + _tmp7;
  const Scalar _tmp20 = _tmp19 + position_vector(0, 0);
  const Scalar _tmp21 = -2 * std::pow(_tmp8, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp3 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp15 + Scalar(0.20999999999999999) * _tmp16;
  const Scalar _tmp24 = _tmp4 * _tmp9;
  const Scalar _tmp25 = _tmp11 * _tmp8;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = -Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = -_tmp23 + _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _tmp29 + position_vector(1, 0);
  const Scalar _tmp31 = Scalar(23.356819799277336) *
                            std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp30), Scalar(2)) +
                        Scalar(3.2276287484906994) *
                            std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp20 - 1), Scalar(2));
  const Scalar _tmp32 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp10 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp34 = -_tmp33;
  const Scalar _tmp35 =
      -Scalar(0.010999999999999999) * _tmp21 - Scalar(0.010999999999999999) * _tmp5;
  const Scalar _tmp36 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp37 = _tmp35 + _tmp36;
  const Scalar _tmp38 = _tmp34 + _tmp37;
  const Scalar _tmp39 = _tmp35 - _tmp36;
  const Scalar _tmp40 = _tmp33 + _tmp39;
  const Scalar _tmp41 = -_tmp22;
  const Scalar _tmp42 = _tmp23 + _tmp27;
  const Scalar _tmp43 = _tmp41 + _tmp42;
  const Scalar _tmp44 = _tmp43 + position_vector(1, 0);
  const Scalar _tmp45 = _tmp14 - _tmp17;
  const Scalar _tmp46 = _tmp45 + _tmp6;
  const Scalar _tmp47 = _tmp46 + position_vector(0, 0);
  const Scalar _tmp48 = Scalar(6.347051629636093) *
                            std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp47), Scalar(2)) +
                        Scalar(70.366961588110229) *
                            std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp44 - 1), Scalar(2));
  const Scalar _tmp49 = _tmp47 + Scalar(-2.5193355532036801);
  const Scalar _tmp50 = Scalar(1.0) / (_tmp49);
  const Scalar _tmp51 = _tmp44 + Scalar(8.3885017487099702);
  const Scalar _tmp52 = _tmp50 * _tmp51;
  const Scalar _tmp53 = _tmp40 * _tmp52;
  const Scalar _tmp54 = _tmp22 + _tmp42;
  const Scalar _tmp55 = _tmp54 + position_vector(1, 0);
  const Scalar _tmp56 = _tmp55 + Scalar(-4.7744369927459998);
  const Scalar _tmp57 = _tmp18 + _tmp6;
  const Scalar _tmp58 = _tmp57 + position_vector(0, 0);
  const Scalar _tmp59 = _tmp58 + Scalar(-2.7171519410699099);
  const Scalar _tmp60 = std::pow(Scalar(std::pow(_tmp56, Scalar(2)) + std::pow(_tmp59, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp61 = _tmp59 * _tmp60;
  const Scalar _tmp62 = _tmp33 + _tmp37;
  const Scalar _tmp63 = _tmp56 * _tmp60;
  const Scalar _tmp64 = -_tmp53 * _tmp61 + _tmp62 * _tmp63;
  const Scalar _tmp65 = Scalar(1.0) / (_tmp52 * _tmp61 - _tmp63);
  const Scalar _tmp66 = _tmp52 * _tmp65;
  const Scalar _tmp67 = _tmp53 + _tmp64 * _tmp66;
  const Scalar _tmp68 = Scalar(1.0) * _tmp43;
  const Scalar _tmp69 = -_tmp68;
  const Scalar _tmp70 = Scalar(1.0) / (_tmp54 + _tmp69);
  const Scalar _tmp71 = Scalar(1.0) * _tmp46;
  const Scalar _tmp72 = _tmp70 * (-_tmp57 + _tmp71);
  const Scalar _tmp73 = _tmp40 * _tmp61 - _tmp61 * _tmp62;
  const Scalar _tmp74 = -_tmp40 + _tmp66 * _tmp73 - _tmp67 * _tmp72;
  const Scalar _tmp75 = _tmp34 + _tmp39;
  const Scalar _tmp76 = _tmp45 + _tmp7;
  const Scalar _tmp77 = _tmp76 + position_vector(0, 0);
  const Scalar _tmp78 = _tmp77 + Scalar(1.9874742031097401);
  const Scalar _tmp79 = _tmp28 + _tmp41;
  const Scalar _tmp80 = _tmp79 + position_vector(1, 0);
  const Scalar _tmp81 = _tmp80 + Scalar(8.3196563720703107);
  const Scalar _tmp82 = std::pow(Scalar(std::pow(_tmp78, Scalar(2)) + std::pow(_tmp81, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp83 = _tmp78 * _tmp82;
  const Scalar _tmp84 = _tmp81 * _tmp82;
  const Scalar _tmp85 = _tmp52 * _tmp83 - _tmp84;
  const Scalar _tmp86 = _tmp65 * _tmp85;
  const Scalar _tmp87 = -_tmp53 * _tmp83 - _tmp64 * _tmp86 + _tmp75 * _tmp84;
  const Scalar _tmp88 = _tmp40 * _tmp83 - _tmp72 * _tmp87 - _tmp73 * _tmp86 - _tmp75 * _tmp83;
  const Scalar _tmp89 = Scalar(1.0) / (_tmp88);
  const Scalar _tmp90 =
      std::sqrt(Scalar(std::pow(_tmp49, Scalar(2)) + std::pow(_tmp51, Scalar(2))));
  const Scalar _tmp91 = Scalar(1.0) / (_tmp90);
  const Scalar _tmp92 = _tmp50 * _tmp90;
  const Scalar _tmp93 = _tmp92 * (-_tmp43 * _tmp49 * _tmp91 + _tmp46 * _tmp51 * _tmp91);
  const Scalar _tmp94 = _tmp54 * _tmp61 - _tmp57 * _tmp63 + _tmp61 * _tmp93;
  const Scalar _tmp95 = -_tmp76 * _tmp84 + _tmp79 * _tmp83 + _tmp83 * _tmp93 - _tmp86 * _tmp94;
  const Scalar _tmp96 = _tmp89 * _tmp95;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp95);
  const Scalar _tmp98 = _tmp88 * _tmp97;
  const Scalar _tmp99 = _tmp98 * (_tmp66 * _tmp94 - _tmp74 * _tmp96 - _tmp93);
  const Scalar _tmp100 = _tmp69 + _tmp79;
  const Scalar _tmp101 = _tmp100 * _tmp72;
  const Scalar _tmp102 = Scalar(1.0) / (-_tmp101 + _tmp71 - _tmp76);
  const Scalar _tmp103 = _tmp100 * _tmp102;
  const Scalar _tmp104 = _tmp74 + _tmp99;
  const Scalar _tmp105 = _tmp87 * _tmp89;
  const Scalar _tmp106 = _tmp103 * _tmp99 - _tmp104 * _tmp105 + _tmp67;
  const Scalar _tmp107 = Scalar(1.0) * _tmp70;
  const Scalar _tmp108 = Scalar(1.0) * _tmp102;
  const Scalar _tmp109 = _tmp20 + Scalar(1.7965602546229);
  const Scalar _tmp110 = _tmp30 + Scalar(-4.83288938413423);
  const Scalar _tmp111 =
      std::pow(Scalar(std::pow(_tmp109, Scalar(2)) + std::pow(_tmp110, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp112 = _tmp109 * _tmp111;
  const Scalar _tmp113 = _tmp112 * fh1;
  const Scalar _tmp114 = _tmp108 * _tmp98;
  const Scalar _tmp115 = Scalar(1.0) * _tmp97;
  const Scalar _tmp116 = _tmp100 * _tmp114 - _tmp115 * _tmp87;
  const Scalar _tmp117 = _tmp110 * _tmp111;
  const Scalar _tmp118 = fh1 * (-_tmp112 * _tmp29 + _tmp117 * _tmp19);
  const Scalar _tmp119 = Scalar(1.0) * _tmp65;
  const Scalar _tmp120 = _tmp119 * _tmp64;
  const Scalar _tmp121 = -_tmp119 * _tmp73 + _tmp120 * _tmp72;
  const Scalar _tmp122 = _tmp98 * (-_tmp119 * _tmp94 - _tmp121 * _tmp96);
  const Scalar _tmp123 = _tmp121 + _tmp122;
  const Scalar _tmp124 = _tmp103 * _tmp122 - _tmp105 * _tmp123 - _tmp120;
  const Scalar _tmp125 = _tmp117 * fh1;
  const Scalar _tmp126 = _tmp100 * _tmp70;
  const Scalar _tmp127 = _tmp38 * fh1;
  const Scalar _tmp128 = _tmp112 * _tmp127 + Scalar(5.1796800000000003) * _tmp13 + _tmp19 * fv1;
  const Scalar _tmp129 = _tmp108 * _tmp72;
  const Scalar _tmp130 = _tmp101 * _tmp108 + Scalar(1.0);
  const Scalar _tmp131 = -_tmp117 * _tmp127 - Scalar(5.1796800000000003) * _tmp26 - _tmp29 * fv1;
  const Scalar _tmp132 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp133 = _tmp68 * _tmp72 + _tmp71;
  const Scalar _tmp134 = _tmp102 * _tmp133;
  const Scalar _tmp135 = 0;
  const Scalar _tmp136 = _tmp135 * _tmp89;
  const Scalar _tmp137 = _tmp70 * (-_tmp100 * _tmp134 - _tmp136 * _tmp87 + _tmp69);
  const Scalar _tmp138 = _tmp61 * _tmp65;
  const Scalar _tmp139 = _tmp85 * _tmp89;
  const Scalar _tmp140 = _tmp135 * _tmp139;
  const Scalar _tmp141 = -_tmp123 * _tmp139 + Scalar(1.0);
  const Scalar _tmp142 = _tmp83 * _tmp89;
  const Scalar _tmp143 = -_tmp104 * _tmp139 - _tmp52;
  const Scalar _tmp144 = _tmp119 * _tmp85 * _tmp97;
  const Scalar _tmp145 = -_tmp113 * _tmp92 * (_tmp104 * _tmp142 + _tmp138 * _tmp143 + Scalar(1.0)) -
                         _tmp118 * _tmp92 * (_tmp115 * _tmp83 - _tmp144 * _tmp61) -
                         _tmp125 * _tmp92 * (_tmp123 * _tmp142 + _tmp138 * _tmp141) -
                         _tmp132 * _tmp92 * (_tmp136 * _tmp83 - _tmp138 * _tmp140);
  const Scalar _tmp146 = Scalar(1.0) / (_tmp145);
  const Scalar _tmp147 =
      std::asinh(_tmp146 * (Scalar(1.0) * _tmp113 * (-_tmp106 * _tmp107 + _tmp108 * _tmp99) +
                            Scalar(1.0) * _tmp118 * (-_tmp107 * _tmp116 + _tmp114) +
                            Scalar(1.0) * _tmp125 * (-_tmp107 * _tmp124 + _tmp108 * _tmp122) +
                            Scalar(1.0) * _tmp128 * (_tmp108 * _tmp126 - _tmp108) +
                            Scalar(1.0) * _tmp131 * (-_tmp107 * _tmp130 + _tmp129) +
                            Scalar(1.0) * _tmp132 *
                                (-_tmp108 * _tmp133 - Scalar(1.0) * _tmp137 + Scalar(1.0))));
  const Scalar _tmp148 = Scalar(9.6622558468725703) * _tmp145;
  const Scalar _tmp149 = _tmp113 * _tmp143 * _tmp65 - _tmp118 * _tmp144 +
                         _tmp125 * _tmp141 * _tmp65 - _tmp132 * _tmp140 * _tmp65;
  const Scalar _tmp150 = Scalar(1.0) / (_tmp149);
  const Scalar _tmp151 = _tmp108 * _tmp128;
  const Scalar _tmp152 =
      std::asinh(_tmp150 * (_tmp106 * _tmp113 * _tmp70 + _tmp116 * _tmp118 * _tmp70 +
                            _tmp124 * _tmp125 * _tmp70 - _tmp126 * _tmp151 +
                            _tmp130 * _tmp131 * _tmp70 + _tmp132 * _tmp137));
  const Scalar _tmp153 = Scalar(22.795248597701466) *
                             std::pow(Scalar(1 - Scalar(0.2094487793051498) * _tmp55), Scalar(2)) +
                         Scalar(7.3829146708599787) *
                             std::pow(Scalar(1 - Scalar(0.36803241838814449) * _tmp58), Scalar(2));
  const Scalar _tmp154 = Scalar(9.6622558468725703) * _tmp149;
  const Scalar _tmp155 = _tmp104 * _tmp113 * _tmp89 + _tmp115 * _tmp118 +
                         _tmp123 * _tmp125 * _tmp89 + _tmp132 * _tmp136;
  const Scalar _tmp156 = Scalar(1.0) / (_tmp155);
  const Scalar _tmp157 =
      std::asinh(_tmp156 * (-_tmp102 * _tmp113 * _tmp99 - _tmp102 * _tmp122 * _tmp125 -
                            _tmp114 * _tmp118 - _tmp129 * _tmp131 + _tmp132 * _tmp134 + _tmp151));
  const Scalar _tmp158 =
      Scalar(3.9500537080266964) *
          std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp77 - 1), Scalar(2)) +
      Scalar(69.216682149330126) *
          std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp80 - 1), Scalar(2));
  const Scalar _tmp159 = Scalar(9.6622558468725703) * _tmp155;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp32 * (-std::sinh(Scalar(1.0) * _tmp1) -
                std::sinh(Scalar(0.1034955) * _tmp0 * (-_tmp1 * _tmp32 - std::sqrt(_tmp31)))) -
      Scalar(8.3701287145205097) *
          std::sqrt(Scalar(Scalar(0.014273672465547484) * _tmp31 +
                           std::pow(Scalar(-Scalar(0.11947247576553975) * _tmp38 -
                                           Scalar(0.11947247576553975) * position_vector(2, 0) + 1),
                                    Scalar(2))));
  _res(1, 0) =
      _tmp148 *
          (-std::sinh(Scalar(1.0) * _tmp147) -
           std::sinh(Scalar(0.1034955) * _tmp146 * (-_tmp147 * _tmp148 - std::sqrt(_tmp48)))) -
      Scalar(8.4690207536791995) *
          std::sqrt(Scalar(Scalar(0.013942273753185566) * _tmp48 +
                           std::pow(Scalar(-Scalar(0.11807740576920533) * _tmp40 -
                                           Scalar(0.11807740576920533) * position_vector(2, 0) + 1),
                                    Scalar(2))));
  _res(2, 0) =
      _tmp154 *
          (-std::sinh(Scalar(1.0) * _tmp152) -
           std::sinh(Scalar(0.1034955) * _tmp150 * (-_tmp152 * _tmp154 - std::sqrt(_tmp153)))) -
      Scalar(8.3641088633305802) *
          std::sqrt(Scalar(Scalar(0.014294226073078949) * _tmp153 +
                           std::pow(Scalar(-Scalar(0.11955846299229073) * _tmp62 -
                                           Scalar(0.11955846299229073) * position_vector(2, 0) + 1),
                                    Scalar(2))));
  _res(3, 0) =
      _tmp159 *
          (-std::sinh(Scalar(1.0) * _tmp157) -
           std::sinh(Scalar(0.1034955) * _tmp156 * (-_tmp157 * _tmp159 - std::sqrt(_tmp158)))) -
      Scalar(8.4718465805053693) *
          std::sqrt(Scalar(Scalar(0.013932974274013004) * _tmp158 +
                           std::pow(Scalar(-Scalar(0.11803802045956635) * _tmp75 -
                                           Scalar(0.11803802045956635) * position_vector(2, 0) + 1),
                                    Scalar(2))));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
