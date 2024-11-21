import { truncate } from "lodash-es";
import { EChart } from "@kbox-labs/react-echarts";
import classes from "./RangeChart.module.css";

/** chart that shows values as points and ranges as bars */
const RangeChart = ({
  title,
  xLabel,
  yLabels,
  values,
  lowerBounds,
  upperBounds,
  tooltip = () => "",
}) => {
  /** difference between upper and lower */
  const upperLowerDiff = lowerBounds.map(
    (_, index) => upperBounds[index] - lowerBounds[index],
  );

  return (
    <EChart
      height={50 + values.length * 40}
      className={classes.chart}
      grid={{ containLabel: true }}
      renderer="svg"
      textStyle={{ fontFamily: "inherit", color: "var(--black)" }}
      title={{
        text: title,
        left: "center",
        textStyle: {
          fontSize: 20,
        },
      }}
      xAxis={{
        type: "value",
        name: xLabel,
        axisLabel: {
          fontSize: 20,
        },
        nameTextStyle: {
          fontSize: 20,
        },
        nameLocation: "middle",
        nameGap: 50,
      }}
      yAxis={{
        type: "category",
        data: yLabels,
        axisLabel: {
          fontSize: 20,
          padding: 20,
          formatter: (value) =>
            truncate(value, {
              length: window.innerWidth / 30,
            }),
        },
        nameTextStyle: {
          fontSize: 20,
        },
      }}
      series={[
        {
          data: lowerBounds,
          type: "bar",
          stack: "group",
          itemStyle: {
            opacity: 0,
          },
        },
        {
          data: upperLowerDiff,
          type: "bar",
          barWidth: 10,
          stack: "group",
          itemStyle: {
            color: "var(--primary)",
            opacity: 1,
          },
          emphasis: {
            itemStyle: {
              color: "color-mix(in srgb, var(--primary), var(--white) 50%)",
            },
          },
        },
        {
          data: values,
          type: "scatter",
          itemStyle: {
            color: "var(--secondary)",
            opacity: 1,
          },
          emphasis: {
            itemStyle: {
              color: "color-mix(in srgb, var(--secondary), var(--white) 50%)",
            },
          },
          symbol: "diamond",
          symbolSize: 20,
        },
      ]}
      tooltip={{
        trigger: "axis",
        axisPointer: {
          type: "shadow",
        },
        formatter: ([series]) =>
          `<div class="${classes.tooltip}">
            ${tooltip(series.dataIndex)}
          </div>`,
      }}
    />
  );
};

export default RangeChart;
