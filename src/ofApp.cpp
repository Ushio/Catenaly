#include "ofApp.h"

#include <functional>
#include <random>


double I(double a, double b, std::function<double(double)> y, std::function<double(double)> dydx) {
	int N = 100;
	double dx = (b - a) / N;
	double sum = 0;
	for (int i = 0; i < N; ++i) {
		double x = a + dx * i;
		double dydx_value = dydx(x);
		sum += y(x) * std::sqrt(1.0 + dydx_value * dydx_value) * dx;
	}
	return sum;
}

double function_y(double x) {
	return 10;
}

double function_dydx(double x) {
	return 10;
}

double function_delta_y(double x) {
	return 0;
}

double function_delta_dydx(double x) {
	return 0;
}


/*
omposite Simpson's rule
n is even
*/
inline double integrate_composite_simpson(std::function<double(double)> f, double a, double b, int n) {
	assert((n & 0x1) == 0);
	double sum = 0;
	double h = (b - a) / n;
	for (int i = 1; i < n; ++i) {
		double c = (i & 0x1) ? 4.0 : 2.0;
		double x_i = a + h * i;
		sum += c * f(x_i);
	}
	sum += f(a) + f(b);
	return sum * h / 3.0;
}

namespace catenary {
	template <int n>
	static double Pow(double v) {
		static_assert(n > 0, "Power can’t be negative");
		double n2 = Pow<n / 2>(v);
		return n2 * n2 * Pow<n & 1>(v);
	}
	template <> double Pow<1>(double v) { return v; }
	template <> double Pow<0>(double v) { return 1.0; }

	// important
	// h: holizontal distance
	// v: vertical distance

	// T = - Sqrt(s^2 - v^2)
	// f(a) = 2*a*sinh(h/(a2)) + T
	// for finding a, f(a) == 0
	inline double f(double a, double h, double sqrt_s2_minus_v2) {
		double T = -sqrt_s2_minus_v2;
		return 2.0 * a * std::sinh(h / (a * 2.0)) + T;
	}
	// f'(a)
	// D[2*a*sinh(h/(2a)) + T, a]
	inline double dfdx(double a, double h) {
		return 2.0 * std::sinh(h / (2.0 * a)) - (h * std::cosh(h / (2.0 * a)) / a);
	}

	// taylor series 2a*sinh(h/(2a)) + T
	inline double f_taylor(double a, double h, double sqrt_s2_minus_v2) {
		double T = -sqrt_s2_minus_v2;
		// T + h + h^3 / (24 a^2) + h^5 / (1920 a^4)
		// T + h * { 1.0 + h^2 / (24 a^2) + h^4 / (1920 a^4) }
		// T + h * { 1.0 + (h^2 / (a^2)) { (1/24 + h^2/(1920 a^2) } }
		// T + h * { 1.0 + (h^2 / (a^2)) { (1/24 + (h^2 / (a^2)) / 1920) } }
		double a2 = a * a;
		double a4 = a2 * a2;
		double h2 = h * h;
		double h2_over_a2 = h2 / a2;
		return std::fma(
			h,
			std::fma(
				h2_over_a2,
				std::fma(h2_over_a2, 1.0 / 1920.0, 1.0 / 24.0),
				1.0),
			T);
	}

	// f(a) = y
	// f^(-1)(y) = a
	inline double f_inverse_taylor(double y, double h, double sqrt_s2_minus_v2) {
		double T = -sqrt_s2_minus_v2;

		// constant
		static const double k_one_over_4_sqrt_15 = 1.0 / (4.0 * std::sqrt(15.0));
		static const double k_sqrt_5 = std::sqrt(5.0);

		// Solve[T+h + h^3/(24 a^2) + h^5/(1920 a^4)==x, a]
		// CForm(sqrt(-(sqrt(5) sqrt(-h^5 (h + 6 T - 6 x)) + 5 h^3)/(h + T - x))/(4 sqrt(15)))
		return std::sqrt(-(k_sqrt_5 * std::sqrt(-(Pow<5>(h) * (h + 6.0 * T - 6.0 * y))) + 5 * Pow<3>(h)) / (h + T - y)) * k_one_over_4_sqrt_15;
	}

	// find a, f(a) = y
	// f^(-1)(y) = a
	inline double f_inverse_newton(double y, double h, double sqrt_s2_minus_v2, int *iteratedCount = nullptr) {
		double T = -sqrt_s2_minus_v2;
		double eps = 1.0e-9;
		double x_cur = f_inverse_taylor(y, h, sqrt_s2_minus_v2);
		double f_value = f(x_cur, h, sqrt_s2_minus_v2) - y;

		if (iteratedCount) {
			*iteratedCount = 1;
		}

		for (int i = 0; i < 100; ++i) {
			if (std::isfinite(x_cur) == false) {
				break;
			}
			if (std::fabs(f_value) < eps) {
				break;
			}
			double move = -f_value / dfdx(x_cur, h);
			x_cur = x_cur + move;
			f_value = f(x_cur, h, sqrt_s2_minus_v2) - y;

			if (iteratedCount) {
				(*iteratedCount)++;
			}
		}
		return x_cur;
	}
	inline double solve_a(double h, double v, double s) {
		return f_inverse_newton(0.0, h, sqrt(s * s - v * v));
	}

	struct Curve {
		double a = 1.0;
		double S = 0.0;
		double T = 0.0;
		double evaluate(double x) const {
			return a * std::cosh((x + S) / a) + T;
		}
		double evaluate_dfdx(double x) const {
			return std::sinh((x + S) / a);
		}
		double x_for_s(double s) const {
			return a * std::asinh(s / a + std::sinh(S / a)) - S;
		}
	};

	/*
	(x1, x2) = (0, 0)
	(x2, y2) = (X, Y)
	length: s Catenary curve
	*/
	inline Curve curve(double X, double Y, double s) {
		Curve c;
		c.a = solve_a(X, std::fabs(Y), s);
		c.S = c.a * std::asinh(Y / (c.a * 2.0 * std::sinh(X / (2.0 * c.a)))) - X * 0.5;
		c.T = -c.a * std::cosh(c.S / c.a);
		return c;
	}

	// double L_actual = c.a * std::sinh((X + c.S) / c.a) - c.a * std::sinh(c.S / c.a);
	// printf("%.3f\n", L_actual);
}

void test() {
	std::mt19937 engine;
	std::uniform_real_distribution<> random_x(0.01, 1000.0);
	std::uniform_real_distribution<> random_y(-1000.0, 1000.0);
	double eps = 0.001;
	for (int i = 0; i < 100000; ++i) {
		double x = random_x(engine);
		double y = random_y(engine);
		double min_s = std::sqrt(x * x + y * y);
		std::uniform_real_distribution<> random_s(min_s + eps, min_s + eps + 1000.0);
		double s = random_s(engine);
		catenary::Curve curve = catenary::curve(x, y, s);

		double numeric_s = integrate_composite_simpson([curve](double x) {
			double dfdx = curve.evaluate_dfdx(x);
			return std::sqrt(1.0 + dfdx * dfdx);
		}, 0, x, 1000);
		
		if (std::fabs(curve.evaluate(0)) > 0.000001) {
			abort();
		}
		if (std::fabs(curve.evaluate(x) - y) > 0.000001) {
			abort();
		}
		if (std::fabs(s - numeric_s) > eps) {
			abort();
		}

		std::uniform_real_distribution<> random_at_s(0.0, s);
		double at_s = random_at_s(engine);
		double at_x = curve.x_for_s(at_s);

		double numeric_at_s = integrate_composite_simpson([curve](double x) {
			double dfdx = curve.evaluate_dfdx(x);
			return std::sqrt(1.0 + dfdx * dfdx);
		}, 0, at_x, 1000);

		if (std::fabs(at_s - numeric_at_s) > eps) {
			abort();
		}
	}
}

//--------------------------------------------------------------
void ofApp::setup() {

	double E = I(0, 10, function_y, function_dydx);

	double E_plus_deltaE = I(0, 10, [](double x) {
		return function_y(x) + function_delta_y(x);
	}, [](double x) {
		return function_dydx(x) + function_delta_dydx(x);
	});

	_camera.setNearClip(0.1f);
	_camera.setFarClip(100.0f);
	_camera.setDistance(5.0f);

	ofSetFrameRate(60);
	ofSetVerticalSync(false);

	_imgui.setup();

	// test();
}
//--------------------------------------------------------------
void ofApp::update() {

}


//--------------------------------------------------------------
void ofApp::draw() {
	float elapsed = ofGetElapsedTimef();

	ofClear(0);
	_camera.begin();
	
	ofPushMatrix();
	ofRotateY(90.0f);
	ofSetColor(64);
	ofDrawGridPlane(1.0f);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(5);
	ofPopMatrix();

	// 元の関数
	//{
	//	ofSetColor(255);

	//	ofPolyline line;
	//	int N = 1000;
	//	for (int i = 1; i < N; ++i) {
	//		float x = ofMap(i, 0, N - 1, 0, 50);
	//		float y = catenary::f(x, 3.0, -0.5);
	//		line.addVertex(x, y);

	//		float dfdx = catenary::dfdx(x, 3.0);
	//		ofVec2f v = ofVec2f(1.0f, dfdx).normalized() * 0.5f;
	//		ofDrawLine(x, y, x + v.x, y + v.y);
	//	}
	//	ofSetColor(255);
	//	line.draw();
	//}

	// テイラー展開ver
	//{
	//	ofPolyline line;
	//	int N = 1000;
	//	for (int i = 1; i < N; ++i) {
	//		float x = ofMap(i, 0, N - 1, 0, 50);
	//		float y = catenary::f_taylor(x, 3.0, -0.5);
	//		line.addVertex(x, y);
	//	}
	//	ofSetColor(255);
	//	line.draw();
	//}

	// 逆関数 テイラー展開ver
	//{
	//	ofPolyline line;
	//	int N = 1000;
	//	for (int i = 0; i < N; ++i) {
	//		float x = ofMap(i, 0, N - 1, 1, 50);
	//		float y = catenary::f_inverse_taylor(x, 3.0, -0.5);
	//		line.addVertex(x, y);
	//	}
	//	ofSetColor(255, 0, 0);
	//	line.draw();
	//}

	// 逆関数 ニュートン
	//{
	//	ofPolyline line;
	//	int N = 1000;
	//	for (int i = 0; i < N; ++i) {
	//		float x = ofMap(i, 0, N - 1, 1, 50);
	//		int iterated = 0;
	//		float y = catenary::f_inverse_newton(x, 3.0, -0.5, &iterated);
	//		line.addVertex(x, y);

	//		ofDrawCircle(x, iterated, 0.1);
	//	}
	//	ofSetColor(255, 0, 255);
	//	line.draw();
	//}

	ofSetColor(128, 255, 128);
	ofDrawSphere(_toX, _toY, 0.05f);

	catenary::Curve curve = catenary::curve(_toX, _toY, _s);
	{
		ofPolyline line;
		int N = 1000;
		for (int i = 0; i < N; ++i) {
			double x = ofMap(i, 0, N - 1, 0, _toX);
			double y = curve.evaluate(x);
			line.addVertex(x, y);
		}
		ofSetColor(255, 255, 0);
		line.draw();
	}

	{
		ofVec3f colors[] = {
			{ 80, 0, 255 },
			{ 0, 204, 128 },
			{ 0, 255, 0 },
			{ 255, 72, 0 },
			{ 224, 0, 0 },
		};

		int N = 10;
		float step = _s / N;
		for (int i = 1; i < N; ++i) {
			auto color = colors[i % (sizeof(colors) / sizeof(colors[i]))];
			ofSetColor(color.x, color.y, color.z);

			double s = step * i;
			double x = curve.x_for_s(s);
			double y = curve.evaluate(x);

			double ofs = 0.6;
			ofLine(x, y, x, y - ofs);
			ofDrawCircle(x, y - ofs, 0.2);
		}
	}

	_camera.end();

	_imgui.begin();
	ImGui::PushStyleColor(ImGuiCol_WindowBg, ofVec4f(0.0f, 0.2f, 0.2f, 1.0f));
	ImGui::SetNextWindowPos(ofVec2f(10, 30), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ofVec2f(500, ofGetHeight() * 0.8), ImGuiSetCond_Once);

	ImGui::Begin("Config Panel");
	
	ImGui::InputFloat("X", &_toX, 0.1f);
	ImGui::InputFloat("Y", &_toY, 0.1f);

	ImGui::InputFloat("s", &_s, 0.1f);

	auto wp = ImGui::GetWindowPos();
	auto ws = ImGui::GetWindowSize();
	ofRectangle win(wp.x, wp.y, ws.x, ws.y);

	ImGui::End();
	ImGui::PopStyleColor();


	_imgui.end();

	if (win.inside(ofGetMouseX(), ofGetMouseY())) {
		_camera.disableMouseInput();
	}
	else {
		_camera.enableMouseInput();
	}
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
