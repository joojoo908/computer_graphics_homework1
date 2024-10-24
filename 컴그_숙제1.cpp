#define _CRT_SECURE_NO_WARNINGS
//#define GLEW_STATIC
//#pragma comment(lib, "glew32s.lib") //왜 이게 되는지는 모름
#pragma comment(lib, "glew32.lib")

#include <gl/glew.h>
#include <gl/freeglut.h>
#include <gl/freeglut_ext.h>
#include <gl/glm/glm.hpp>
#include <gl/glm/ext.hpp>
#include <gl/glm/gtc/matrix_transform.hpp>

#include<iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <stdlib.h>
#include <stdio.h>

struct Color {
	float r, g, b;
};
struct Vertex {
	float x, y, z;
};
struct Face {
	unsigned int v1, v2, v3;
};
struct Model {
	std::vector<Vertex> vertices;
	std::vector<Face> faces;
};

void make_vertexShaders();
void make_fragmentShaders();
GLuint make_shaderProgram();
void InitBuffer(const Model& model);
bool loadOBJ(const std::string& filename, Model& model);
GLvoid drawScene();
void Keyboard(unsigned char key, int x, int y);
void SpecialKeyboard(int key, int x, int y);
void Mouse(int button, int state, int x, int y);
void Motion(int x, int y);
void TimerFunction(int value);
GLvoid Reshape(int w, int h);
char* filetobuf(const char* file);

void clear_mat();
void move(glm::mat4& trans, Vertex move);
void spin(glm::mat4& trans, Vertex spin);
void vector_spin(std::vector<Model>& m, int num, Vertex mid, Vertex spin);
void mid_spin(glm::mat4& trans, Vertex mid, Vertex spin);
void scale(glm::mat4& trans);

Vertex mousevec(int x, int y);
float sign(Vertex v1, Vertex v2, Vertex v3);
bool dot_in_tri(Vertex p, Vertex v1, Vertex v2, Vertex v3);

Vertex middle_Vertex(std::vector<Vertex> vertices) {
	Vertex m_v = { 0,0,0 };
	for (int i = 0; i < vertices.size(); i++) {
		m_v.x += vertices[i].x;
		m_v.y += vertices[i].y;
		m_v.z += vertices[i].z;
	}
	m_v.x /= vertices.size();
	m_v.y /= vertices.size();
	m_v.z /= vertices.size();
	return m_v;
}


void change_line(glm::mat4 t , Vertex &p1, Vertex &p2) {
	glm::mat4 inv = glm::inverse(t);
	glm::vec4 vertexPosition(p1.x, p1.y, p1.z, 1.0f);
	glm::vec4 vertexPosition2(p2.x, p2.y, p2.z, 1.0f);

	glm::vec4 transformedPosition = inv * vertexPosition;
	glm::vec4 transformedPosition2 = inv * vertexPosition2;

	p1.x = transformedPosition.x;
	p1.y = transformedPosition.y;
	p1.z = transformedPosition.z;

	p2.x = transformedPosition2.x;
	p2.y = transformedPosition2.y;
	p2.z = transformedPosition2.z;
}
bool cut_check(Model m, glm::mat4 t, Vertex &p1, Vertex &p2) {
	glm::mat4 inv = glm::inverse(t);
	glm::vec4 vertexPosition(p1.x, p1.y, p1.z, 1.0f);
	glm::vec4 vertexPosition2(p2.x, p2.y, p2.z, 1.0f);
	glm::vec4 transformedPosition = inv * vertexPosition;
	glm::vec4 transformedPosition2 = inv * vertexPosition2;

	p1.x = transformedPosition.x;
	p1.y = transformedPosition.y;
	p1.z = transformedPosition.z;

	p2.x = transformedPosition2.x;
	p2.y = transformedPosition2.y;
	p2.z = transformedPosition2.z;

	Vertex m_v = middle_Vertex(m.vertices);

	return 0;
}
bool ck_line(Vertex AP1, Vertex AP2,Vertex BP1, Vertex BP2 , Vertex &ip){
	float t;
	float s;
	float under = (BP2.y - BP1.y) * (AP2.x - AP1.x) - (BP2.x - BP1.x) * (AP2.y - AP1.y);
	if (under == 0) return false;

	float _t = (BP2.x - BP1.x) * (AP1.y - BP1.y) - (BP2.y - BP1.y) * (AP1.x - BP1.x);
	float _s = (AP2.x - AP1.x) * (AP1.y - BP1.y) - (AP2.y - AP1.y) * (AP1.x - BP1.x);

	t = _t / under;
	s = _s / under;

	if (t < 0.0 || t>1.0 || s < 0.0 || s>1.0) return false;
	if (_t == 0 && _s == 0) return false;

	ip.x = BP1.x + s * (BP2.x - BP1.x);
	ip.y = BP1.y + s * (BP2.y - BP1.y);
	ip.z = BP1.z + s * (BP2.z - BP1.z);

	return true;
}

std::vector < Vertex > dots; //디버깅용

bool model_cut(std::vector<Model> &m, int modelnum , glm::mat4 t, Vertex p1 ,Vertex p2) {
	change_line(t, p1, p2);

	dots.clear();
	Vertex ck = {0,0,0};

	Model up;
	Model down;

	for (int i = 0; i < m[modelnum].faces.size(); i++) {

		if (ck_line(p1, p2, m[modelnum].vertices[m[modelnum].faces[i].v1], m[modelnum].vertices[m[modelnum].faces[i].v2], ck ) ) {
			//std::cout << m[modelnum].faces[i].v1 << " , " << m[modelnum].faces[i].v2 << std::endl;
			dots.push_back(ck);
		}
		if (
			ck_line(p1, p2, m[modelnum].vertices[m[modelnum].faces[i].v2], m[modelnum].vertices[m[modelnum].faces[i].v3], ck)) {
			//std::cout << m[modelnum].faces[i].v2 << " , " << m[modelnum].faces[i].v3 << std::endl;
			dots.push_back(ck);
		}
		if (
			ck_line(p1, p2, m[modelnum].vertices[m[modelnum].faces[i].v3], m[modelnum].vertices[m[modelnum].faces[i].v1], ck)) {
			//std::cout << m[modelnum].faces[i].v3 << " , " << m[modelnum].faces[i].v1 << std::endl;
			dots.push_back(ck);
		}

	}
	for (int i = 0; i < dots.size(); i++) {
		std::cout << dots[i].x <<", " << dots[i].y <<"," << dots[i].z << std::endl;
	}


	//m[modelnum].vertices[0] = p2;
	return 0;
}
//나중에 지울지도
void dots_spin(std::vector < Vertex > &vertices, Vertex mid, Vertex spin) {
	glm::mat4 Rx = glm::mat4(1.0f); //--- 이동 행렬 선언
	glm::mat4 Ry = glm::mat4(1.0f); //--- 회전 행렬 선언
	glm::mat4 Rz = glm::mat4(1.0f); //--- 회전 행렬 선언

	Rx = glm::rotate(Rx, glm::radians(spin.x), glm::vec3(1.0, 0.0, 0.0));
	Ry = glm::rotate(Ry, glm::radians(spin.y), glm::vec3(0.0, 1.0, 0.0));
	Rz = glm::rotate(Rz, glm::radians(spin.z), glm::vec3(0.0, 1.0, 0.0));

	glm::mat4 transformMatrix(1.0f);
	transformMatrix = glm::translate(transformMatrix, glm::vec3(-mid.x, -mid.y, -mid.z));
	glm::mat4 inv = glm::inverse(transformMatrix);

	for (size_t i = 0; i <vertices.size(); i++) {
		glm::vec4 vertexPosition(vertices[i].x, vertices[i].y, vertices[i].z, 1.0f);
		glm::vec4 transformedPosition = inv * Rx * Ry * Rz * transformMatrix * vertexPosition;

		// 변환된 위치를 다시 model에 저장
		vertices[i].x = transformedPosition.x;
		vertices[i].y = transformedPosition.y;
		vertices[i].z = transformedPosition.z;
		//std::cout << spin.x << std::endl;
	}
}

//--- 필요한 변수 선언-----------------------------------------------------------------------------------
GLint width = 1000, height = 800;
GLuint shaderProgramID; //--- 세이더 프로그램 이름
GLuint vertexShader; //--- 버텍스 세이더 객체
GLuint fragmentShader; //--- 프래그먼트 세이더 객체
GLint result;
GLchar errorLog[512];
GLuint shaderID;
GLuint vao, vbo[2], ebo;
Vertex mouse;
//---------------------------------------------------------------------------------------------------

Model mid_line = { {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1} }, };
Model slice;

std::vector<Model> models;

Model model;
glm::mat4 trans(1.0f);
Vertex sp = { 30,-30,0 };
Vertex mid = { -0.5,-0.5,-0.5 };
Vertex v0 = middle_Vertex(model.vertices);
Model dia;
std::vector<Color> colors = {
	{1.0, 0.0, 0.0}, // 빨간색
	{0.0, 1.0, 0.0}, // 초록색
	{0.0, 0.0, 1.0}, // 파란색
	{1.0, 1.0, 0.0}, // 노란색
	{0.0, 1.0, 1.0}, // 시안
	{1.0, 0.0, 1.0}, // 마젠타
	{0.0, 0.0, 0.0}, // 검은색
	{0.5, 0.5, 0.5},  // 흰색
	////////
	{ 1.0, 0.0, 0.0 }, // 빨간색
	{0.0, 1.0, 0.0}, // 초록색
	{0.0, 0.0, 1.0}, // 파란색
	{1.0, 1.0, 0.0}, // 노란색
	{0.0, 1.0, 1.0}, // 시안
	{1.0, 0.0, 1.0} // 마젠타
};
bool cmdc = 1;
bool cmdh = 1;
bool cmdw;
int spin2[2] = { 0,0 };
Vertex spin_model;
Vertex move_model;

//--------------------- 메인 함수----------------------------------------------------------------------------
int main(int argc, char** argv) //--- 윈도우 출력하고 콜백함수 설정
{
	//--- 윈도우 생성하기
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(300, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow("Example1"); //--- GLEW 초기화하기
	glewExperimental = GL_TRUE;
	glewInit();
	//--- 세이더 읽어와서 세이더 프로그램 만들기

	loadOBJ("cube.obj", model);
	//loadOBJ("dia.obj", dia);
	
	models.push_back(model);
	//std::cout << "Vertices in dia: " << model.vertices.size() << std::endl;
	//std::cout << "Faces in dia: " << model.faces.size() << std::endl;

	/*for (int i = 0; i < model.faces.size(); i++) {
		std::cout << model.faces[i].v1 << "  " << model.faces[i].v2 << " " << model.faces[i].v3 << std::endl;
	}*/

	Vertex v0 =middle_Vertex(model.vertices);
	scale(trans);
	move(trans, mid);
	//mid_spin(trans, v0, sp);

	//InitBuffer(model);
	shaderProgramID = make_shaderProgram();
	//--- 세이더 프로그램 만들기
	glutDisplayFunc(drawScene); //--- 출력 콜백 함수
	glutKeyboardFunc(Keyboard);
	glutSpecialFunc(SpecialKeyboard);
	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);
	glutTimerFunc(100, TimerFunction, 1);
	glutReshapeFunc(Reshape);
	glutMainLoop();
}

GLvoid drawScene() //--- 콜백 함수: 그리기 콜백 함수
{
	GLfloat rColor, gColor, bColor;
	rColor = gColor = bColor = 1.0;
	glClearColor(rColor, gColor, bColor, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	//glDisable(GL_CULL_FACE);

	clear_mat();
	/*InitBuffer(mid_line);
	glPointSize(5.0);
	glDrawArrays(GL_LINES, 0, 4);*/
	if (slice.vertices.size() == 2) {
		InitBuffer(slice);
		glPointSize(5.0);
		glDrawArrays(GL_LINES, 0, 2);
		//glm::mat4 inv = glm::inverse(trans); // 역행렬 계산
		//unsigned int transformLocation2 = glGetUniformLocation(shaderID, "transform");
		//glUniformMatrix4fv(transformLocation2, 1, GL_FALSE, glm::value_ptr(inv));
		//InitBuffer(slice);
		//glPointSize(5.0);
		//glDrawArrays(GL_LINES, 0, 2);
	}


	Model cpy;
	for (int i = 0; i < models.size(); i++) {
		cpy = models[i];
		unsigned int transformLocation = glGetUniformLocation(shaderID, "transform");
		glUniformMatrix4fv(transformLocation, 1, GL_FALSE, glm::value_ptr(trans));
		InitBuffer(cpy);
		if (cmdw) {
			glDrawElements(GL_TRIANGLES, cpy.faces.size() * 3, GL_UNSIGNED_INT, 0);
		}
		else {
			glDrawElements(GL_LINE_STRIP, cpy.faces.size() * 3, GL_UNSIGNED_INT, 0);
		}
	}

	//glm::mat4 inv = glm::inverse(trans); // 역행렬 계산
	//clear_mat();
	Model mod = { dots };
	InitBuffer(mod);
	glPointSize(5.0);
	glDrawArrays(GL_POINTS, 0, dots.size());

	glutSwapBuffers(); // 화면에 출력하기
	//glutSwapBuffers(); // 화면에 출력하기
}

void Keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'c':
		cmdc = 1;
		break;
	case 'p':
		cmdc = 0;
		break;
	case 'h':
		cmdh = !cmdh;
		break;
	case 'w':
		cmdw = !cmdw;
		break;
	case 'x':
		spin_model.x = (spin_model.x == 0) ? 1 : 0;
		break;
	case 'X':
		spin_model.x = (spin_model.x == 0) ? -1 : 0;
		break;
	case 'y':
		spin_model.y = (spin_model.y == 0) ? 1 : 0;
		break;
	case 'Y':
		spin2[1] = (spin2[1] == 0) ? -1 : 0;
		break;
	case 's':
		move_model.x = 0;
		move_model.y = 0;
		spin2[0] = 0;
		spin2[1] = 0;
		break;
	}
}
void SpecialKeyboard(int key, int x, int y)
{
	switch (key) {
	case GLUT_KEY_UP:
		move_model = {0,0.1,0};
		move(trans, move_model);
		break;
	case GLUT_KEY_DOWN:
		move_model = { 0,-0.1,0 };
		move(trans, move_model);
		break;
	case GLUT_KEY_RIGHT:
		move_model = {0.1,0,0 };
		move(trans, move_model);
		break;
	case GLUT_KEY_LEFT:
		move_model = { -0.1,0,0 };
		move(trans, move_model);
		break;
	}

}

void Mouse(int button, int state, int x, int y)
{
	
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		mouse = mousevec(x, y);
		slice.vertices.push_back(mouse);
	}
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
		//mouse = mousevec(x, y);
		model_cut(models, 0, trans, slice.vertices[0], slice.vertices[1]);
		slice.vertices.clear();
		//change_line(trans, slice.vertices[0], slice.vertices[1]);
	}
	std::cout << "x = " << mouse.x << " y = " << mouse.y << std::endl;
}
void Motion(int x, int y)
{
	mouse = mousevec(x, y);
	if (slice.vertices.size() == 1) {
		slice.vertices.push_back(mouse);
	}
	else if (slice.vertices.size() == 2) {
		slice.vertices[1] = mouse;
		//change_line(trans, slice.vertices[0], slice.vertices[1]);
		//std::cout << slice.vertices[0].x << std::endl;
	}

}

void TimerFunction(int value)
{
	//mid_spin(trans, v0, spin_model);

	Vertex v0 = middle_Vertex(model.vertices);
	//mid_spin(trans, v0, spin_model);
	vector_spin(models, 0, v0, spin_model);
	dots_spin(dots, v0, spin_model);


	//spin_model.x += spin2[0];
	//spin_model.y += spin2[1];
	glutPostRedisplay(); // 화면 재 출력
	glutTimerFunc(10, TimerFunction, 1); // 타이머함수 재 설정
}


//-------------------------------------------------------------------------
bool loadOBJ(const std::string& filename, Model& model) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "파일을 열 수 없습니다: " << filename << std::endl;
		return false;
	}

	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string prefix;
		iss >> prefix;

		if (prefix == "v") {
			Vertex vertex;
			iss >> vertex.x >> vertex.y >> vertex.z;
			model.vertices.push_back(vertex);
		}
		else if (prefix == "vn") {

		}
		else if (prefix == "f") {
			Face face;
			unsigned int v[3], n[3];
			char slash;

			for (int i = 0; i < 3; ++i) {
				if (iss >> v[i]) {
					// v/n 형식일 때 (v/n이 있을 경우)
					if (iss.peek() == '/') {
						iss >> slash; // 첫 번째 '/' 읽기
						if (iss.peek() == '/') {
							iss >> slash; // 두 번째 '/' 읽기 (법선 인덱스가 없을 때)
							iss >> n[i];
						}
					}
				}
			}
			face.v1 = v[0] - 1;
			face.v2 = v[1] - 1;
			face.v3 = v[2] - 1;
			model.faces.push_back(face);
		}
	}

	return true;
}

void InitBuffer(const Model& model) {
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// VBO 생성 및 정점 데이터 전송
	glGenBuffers(1, &vbo[0]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, model.vertices.size() * sizeof(Vertex), model.vertices.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)0);
	glEnableVertexAttribArray(0);

	// 색상 VBO 생성 및 색상 데이터 전송
	glGenBuffers(1, &vbo[1]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(Color), colors.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Color), (GLvoid*)0);
	glEnableVertexAttribArray(1);

	// EBO 생성 및 데이터 전송
	if (!model.faces.empty()) {
		glGenBuffers(1, &ebo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
		std::vector<unsigned int> indices;
		for (const auto& face : model.faces) {
			indices.push_back(face.v1);
			indices.push_back(face.v2);
			indices.push_back(face.v3);
		}
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
	}

	//glBindVertexArray(0);
	glBindVertexArray(vao);
}

void clear_mat() {
	glm::mat4 transformMatrix(1.0f);
	unsigned int transformLocation = glGetUniformLocation(shaderID, "transform");
	glUniformMatrix4fv(transformLocation, 1, GL_FALSE, glm::value_ptr(transformMatrix));
}
void move(glm::mat4& trans, Vertex move) {
	glm::mat4 transformMatrix(1.0f);
	transformMatrix = glm::translate(transformMatrix, glm::vec3(move.x, move.y, move.z));

	trans = trans * transformMatrix;

}
void spin(glm::mat4& trans, Vertex spin) {

	glm::mat4 Rx = glm::mat4(1.0f); //--- 이동 행렬 선언
	glm::mat4 Ry = glm::mat4(1.0f); //--- 회전 행렬 선언
	glm::mat4 Rz = glm::mat4(1.0f); //--- 회전 행렬 선언

	Rx = glm::rotate(Rx, glm::radians(spin.x), glm::vec3(1.0, 0.0, 0.0));
	Ry = glm::rotate(Ry, glm::radians(spin.y), glm::vec3(0.0, 1.0, 0.0));
	Rz = glm::rotate(Rz, glm::radians(spin.z), glm::vec3(0.0, 1.0, 0.0));

	trans *= Rx;
	trans *= Ry;
	trans *= Rz;

}
void mid_spin(glm::mat4& trans, Vertex mid, Vertex spin) {
	glm::mat4 transformMatrix(1.0f);
	transformMatrix = glm::translate(transformMatrix, glm::vec3(-mid.x, -mid.y, -mid.z));
	glm::mat4 inv = glm::inverse(transformMatrix);

	glm::mat4 Rx = glm::mat4(1.0f); //--- 이동 행렬 선언
	glm::mat4 Ry = glm::mat4(1.0f); //--- 회전 행렬 선언
	glm::mat4 Rz = glm::mat4(1.0f); //--- 회전 행렬 선언

	Rx = glm::rotate(Rx, glm::radians(spin.x), glm::vec3(1.0, 0.0, 0.0));
	Ry = glm::rotate(Ry, glm::radians(spin.y), glm::vec3(0.0, 1.0, 0.0));
	Rz = glm::rotate(Rz, glm::radians(spin.z), glm::vec3(0.0, 1.0, 0.0));


	trans = trans * inv *Rx*Ry*Rz *transformMatrix;
}
void vector_spin(std::vector<Model>& m, int num, Vertex mid,Vertex spin) {
	glm::mat4 Rx = glm::mat4(1.0f); //--- 이동 행렬 선언
	glm::mat4 Ry = glm::mat4(1.0f); //--- 회전 행렬 선언
	glm::mat4 Rz = glm::mat4(1.0f); //--- 회전 행렬 선언

	Rx = glm::rotate(Rx, glm::radians(spin.x), glm::vec3(1.0, 0.0, 0.0));
	Ry = glm::rotate(Ry, glm::radians(spin.y), glm::vec3(0.0, 1.0, 0.0));
	Rz = glm::rotate(Rz, glm::radians(spin.z), glm::vec3(0.0, 1.0, 0.0));

	glm::mat4 transformMatrix(1.0f);
	transformMatrix = glm::translate(transformMatrix, glm::vec3(-mid.x, -mid.y, -mid.z));
	glm::mat4 inv = glm::inverse(transformMatrix);

	for (size_t i = 0; i < m[num].vertices.size(); i++) {
		glm::vec4 vertexPosition(m[num].vertices[i].x, m[num].vertices[i].y, m[num].vertices[i].z, 1.0f);
		glm::vec4 transformedPosition = inv * Rx * Ry * Rz * transformMatrix * vertexPosition;

		// 변환된 위치를 다시 model에 저장
		m[num].vertices[i].x = transformedPosition.x;
		m[num].vertices[i].y = transformedPosition.y;
		m[num].vertices[i].z = transformedPosition.z;
		//std::cout << spin.x << std::endl;
	}
}
void scale(glm::mat4& trans) {
	glm::mat4 sc = glm::mat4(1.0f);
	sc = glm::scale(sc, glm::vec3(0.5, 0.5, 0.5));
	trans *= sc;
}

Vertex mousevec(int x, int y) {
	Vertex alpa = { ( ((float)x-(width/2))/ (width / 2)) ,  (height/2 -(float)y)/(height / 2)  , 0};
	return alpa;
}
float sign(Vertex v1, Vertex v2, Vertex v3) {
	return (v1.x - v3.x) * (v2.y - v3.y) - (v2.x - v3.y) * (v1.y - v3.y);
}
bool dot_in_tri(Vertex p, Vertex v1, Vertex v2, Vertex v3) {
	float d1, d2, d3;
	bool neg, pos;

	d1 = sign(p, v1, v2);
	d2 = sign(p, v2, v3);
	d3 = sign(p, v3, v1);

	neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
	pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

	return !(neg && pos);
}

//언젠간 쓰겄지
/*void spin(glm::mat4 &trans, Model &model) {
	// 회전 변환 행렬 생성 (예: z축을 기준으로 45도 회전)
	float angle = glm::radians(45.0f); // 45도를 라디안으로 변환
	glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0.0f, 0.0f, 1.0f));

	// 각 버텍스에 대해 회전 변환 적용
	for (size_t i = 0; i < model.vertices.size(); ++i) {
		glm::vec4 vertexPosition(model.vertices[i].x, model.vertices[i].y, model.vertices[i].z, 1.0f);
		glm::vec4 transformedPosition = rotationMatrix * vertexPosition;

		// 변환된 위치를 다시 model에 저장
		model.vertices[i].x = transformedPosition.x;
		model.vertices[i].y = transformedPosition.y;
		model.vertices[i].z = transformedPosition.z;
	}
}
*/
//-------------------------------------------------------------------------------
//--- 다시그리기 콜백 함수
GLvoid Reshape(int w, int h) //--- 콜백 함수: 다시 그리기 콜백 함수
{
	glViewport(0, 0, w, h);
}
char* filetobuf(const char* file)
{
	FILE* fptr;
	long length;
	char* buf;
	fptr = fopen(file, "rb"); // Open file for reading
	if (!fptr) // Return NULL on failure
		return NULL;
	fseek(fptr, 0, SEEK_END); // Seek to the end of the file
	length = ftell(fptr); // Find out how many bytes into the file we are
	buf = (char*)malloc(length + 1); // Allocate a buffer for the entire length of the file and a null terminator
	fseek(fptr, 0, SEEK_SET); // Go back to the beginning of the file
	fread(buf, length, 1, fptr); // Read the contents of the file in to the buffer
	fclose(fptr); // Close the file
	buf[length] = 0; // Null terminator
	return buf; // Return the buffer
}
void make_vertexShaders()
{
	GLchar* vertexSource;
	//--- 버텍스 세이더 읽어 저장하고 컴파일 하기
	//--- filetobuf: 사용자정의 함수로 텍스트를 읽어서 문자열에 저장하는 함수
	vertexSource = filetobuf("vertex.glsl");
	vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexSource, NULL);
	glCompileShader(vertexShader);
	GLint result;
	GLchar errorLog[512];
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &result);
	if (!result)
	{
		glGetShaderInfoLog(vertexShader, 512, NULL, errorLog);
		std::cerr << "ERROR: vertex shader 컴파일 실패\n" << errorLog << std::endl;
		return;
	}
}
void make_fragmentShaders()
{
	GLchar* fragmentSource;
	//--- 프래그먼트 세이더 읽어 저장하고 컴파일하기
	fragmentSource = filetobuf("fragment.glsl"); // 프래그세이더 읽어오기
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
	glCompileShader(fragmentShader);
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &result);
	if (!result)
	{
		glGetShaderInfoLog(fragmentShader, 512, NULL, errorLog);
		std::cerr << "ERROR: frag_shader 컴파일 실패\n" << errorLog << std::endl;
		return;
	}
}
GLuint make_shaderProgram()
{
	make_vertexShaders(); //--- 버텍스 세이더 만들기
	make_fragmentShaders(); //--- 프래그먼트 세이더 만들기

	shaderID = glCreateProgram(); //--- 세이더 프로그램 만들기
	glAttachShader(shaderID, vertexShader); //--- 세이더 프로그램에 버텍스 세이더 붙이기
	glAttachShader(shaderID, fragmentShader); //--- 세이더 프로그램에 프래그먼트 세이더 붙이기
	glLinkProgram(shaderID); //--- 세이더 프로그램 링크하기
	glDeleteShader(vertexShader); //--- 세이더 객체를 세이더 프로그램에 링크했음으로, 세이더 객체 자체는 삭제 가능
	glDeleteShader(fragmentShader);
	glGetProgramiv(shaderID, GL_LINK_STATUS, &result); // ---세이더가 잘 연결되었는지 체크하기
	if (!result) {
		glGetProgramInfoLog(shaderID, 512, NULL, errorLog);
		std::cerr << "ERROR: shader program 연결 실패\n" << errorLog << std::endl;
		return false;
	}
	glUseProgram(shaderID); //--- 만들어진 세이더 프로그램 사용하기
	//--- 여러 개의 세이더프로그램 만들 수 있고, 그 중 한개의 프로그램을 사용하려면
	//--- glUseProgram 함수를 호출하여 사용 할 특정 프로그램을 지정한다.
	//--- 사용하기 직전에 호출할 수 있다.
	return shaderID;
}