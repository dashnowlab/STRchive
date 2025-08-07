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
          /** https://echarts.apache.org/en/option.html#series-boxplot.data */
          data: values.map((value, index) => [
            lowerBounds[index],
            value,
            value,
            value,
            upperBounds[index],
          ]),
          type: "boxplot",
          itemStyle: {
            borderWidth: 3,
            borderColor: "var(--primary)",
          },
        },
        {
          data: values,
          type: "scatter",
          itemStyle: {
            color: "var(--black)",
            opacity: 1,
          },
          emphasis: {
            itemStyle: {
              color: "var(--black)",
            },
          },
          symbolSize: 12,
        },
      ]}
      tooltip={{
        trigger: "axis",
        formatter: ([series]) =>
          `<div class="col ${classes.tooltip}">
            ${tooltip(series.dataIndex)}
          </div>`,
      }}
    />
  );
};

export default RangeChart;
