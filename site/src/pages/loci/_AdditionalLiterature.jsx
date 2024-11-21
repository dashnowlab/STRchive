import ShowMoreItems from "@/components/ShowMoreItems";
import Citation from "./_Citation";

/** additional literature section */
const AdditionalLiterature = ({ additionalLiterature }) => (
  <ShowMoreItems
    items={additionalLiterature.map((citation, index) => (
      <Citation key={index} {...citation} />
    ))}
  />
);

export default AdditionalLiterature;
