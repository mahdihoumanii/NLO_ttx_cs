#pragma once

namespace ttbar_qtsub {

double color_correlation_qq_1_2(double shat, double mt, double mu, double gs, double c);
double color_correlation_qq_1_3(double shat, double mt, double mu, double gs, double c);
double color_correlation_qq_1_4(double shat, double mt, double mu, double gs, double c);
double color_correlation_qq_2_3(double shat, double mt, double mu, double gs, double c);
double color_correlation_qq_2_4(double shat, double mt, double mu, double gs, double c);
double color_correlation_qq_3_4(double shat, double mt, double mu, double gs, double c);

double color_correlation_gg_1_2(double shat, double mt, double mu, double gs, double c);
double color_correlation_gg_1_3(double shat, double mt, double mu, double gs, double c);
double color_correlation_gg_1_4(double shat, double mt, double mu, double gs, double c);
double color_correlation_gg_2_3(double shat, double mt, double mu, double gs, double c);
double color_correlation_gg_2_4(double shat, double mt, double mu, double gs, double c);
double color_correlation_gg_3_4(double shat, double mt, double mu, double gs, double c);

// Dispatchers: handle diagonal (i==j) and symmetric off-diagonal pairs.
// QQ diagonal: CF * born_qq  (all four external legs are quarks)
// GG diagonal: emitter Casimir times born_gg
//   legs 1,2 -> CA * born_gg
//   legs 3,4 -> CF * born_gg
// Throws std::invalid_argument for indices outside 1..4.
double color_correlation_qq(double shat, double mt, double mu, double gs, double c, int i, int j);
double color_correlation_gg(double shat, double mt, double mu, double gs, double c, int i, int j);

}  // namespace ttbar_qtsub
