require "./spec_helper"

describe Glushkov do
  # TODO: Write tests

  it "works" do
    a = Matrix(Int32).by_rows([[1, 3, 1], [1, 0, 0]])
    b = Matrix(Int32).by_rows([[0, 0, 5], [7, 5, 0]])
    puts (a + b).inspect

    a = Matrix(Int32).by_rows([[2, 3, 4], [1, 0, 0]])
    b = Matrix(Int32).by_rows([[0, 1000], [1, 100], [0, 10]])
    puts (a * b).inspect
  end

  it "det 4" do
    a = Matrix(Int32).by_rows([
      [5, -7, 2, 2],
      [0, 3, 0, -4],
      [-5, -8, 0, 3],
      [0, 5, 0, -6],
    ])
    a.det.should eq(20)
  end

  it "det 3" do
    a = Matrix(Int32).by_rows([
      [-7, 5, -5],
      [0, -1, -5],
      [0, -3, 8],
    ])
    a.det.should eq(161)
  end

  it "det 2" do
    a = Matrix(Int32).by_rows([
      [1, 2],
      [3, 4],
    ])
    a.det.should eq(-2)
  end

  it "==" do
    a = Matrix(Int32).by_rows([
      [1, 2],
      [3, 4],
    ])
    b = Matrix(Int32).by_rows([
      [1, 2],
      [3, 4],
    ])
    (a == b).should be_true
  end

  it "sub" do
    a = Matrix(Int32).by_rows([
      [1, 4, 7],
      [3, 0, 5],
      [-1, 9, 11],
    ])

    a.sub(1, 2).should eq(Matrix(Int32).by_rows([
      [1, 4],
      [-1, 9],
    ]))
  end

  it "m" do
    a = Matrix(Int32).by_rows([
      [1, 4, 7],
      [3, 0, 5],
      [-1, 9, 11],
    ])
    a.m(1, 2).should eq(13)
  end

  it "cofactor" do
    a = Matrix(Int32).by_rows([
      [1, 4, 7],
      [3, 0, 5],
      [-1, 9, 11],
    ])
    a.cofactor(1, 2).should eq(-13)
  end
end
